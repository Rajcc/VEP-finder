import requests 
import re
from flask import Flask
from flask import render_template,request
import time
from Bio.Blast import NCBIXML
from io import StringIO
import io
import mutations
from mutations import detect_mutation


app = Flask(__name__)
# whenever user opens root folder gene.html page will show
@app.route('/')
def index():
    return render_template('gene.html')

# whenever user fills the form and sends or post data it is send to gene file which is this current flask file and whenever 
@app.route('/gene',methods=['POST'])

def get_gene():
    gene_id=request.form.get('enter')
    drop=request.form.get('cas')
    print("Dropdown value selected:", drop)
    sequence=None
    vep_results = []
    clinvar_results = []


    if gene_id:
  # Get gene info
        info_url = f"https://rest.ensembl.org/lookup/symbol/{drop}/{gene_id}?expand=1"##look up id is is that can be manipilated by user
        info_res = requests.get(info_url, headers={"Content-Type": "application/json"})

        if not info_res.ok:
            return render_template("error.html", message="❌ Gene not found for selected species.")
    
        if info_res.ok:
            gene_info = info_res.json() ##here we are creating a dictionary of data in json
            strand = gene_info.get("strand") ##here the data in json is in "strand" which we are storing in strand variable
            start = gene_info.get("start")
            end = gene_info.get("end")
            chromosome = gene_info.get("seq_region_name")
            ensembl_id=gene_info.get("id")
     
            server = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/x-fasta"
            
            r = requests.get(server, headers={ "Content-Type" : "text/plain"})
            if r.ok:
                sequence=r.text
            else:
               return render_template("error.html", message="Could not retrieve fasta sequence")
            
        else:
            print("please enter valid ensembl gene id")

        region=f"{chromosome}:{start}-{end}"  
        overlap_url=f"https://rest.ensembl.org/overlap/region/{drop}/{region}?feature=variation"
        print(overlap_url)
        overlap_res=requests.get(overlap_url,headers={"Content-Type":"application/json"})

        if overlap_res.ok:
            variants=overlap_res.json()
            for variant in variants[:10]:#1o for performance
                var_id=variant.get("id")
                pas=variant.get("start")
                vep_url=f"https://rest.ensembl.org/vep/{drop}/id/{var_id}?"
                print(vep_url)
                vep_res=requests.get(vep_url,headers={"content-type":"application/json"})

                if vep_res.ok:
                    vep=vep_res.json()
                    if vep:
                        consequence=vep[0].get("most_severe_consequence","unknown")
                        transcript=vep[0].get("transcript_consequences",[{}]) [0]
                        gene_symbol=transcript.get("gene_symbol","")
                        impact=transcript.get("impact","")
                        vep_results.append({
                            "id": var_id,
                            "position": pas,
                            "consequence": consequence,
                            "impact":impact,
                            "gene": gene_symbol
                        })
                        
                         # Step 2: ClinVar Lookup (only if variant has rsID)
                        if var_id.startswith("rs"):
                            clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/v0/beta/refsnp{var_id[2:]}" #api for clinvar 
                            ##results for snp
                            clinvar_res=requests.get(clinvar_url)

                            if clinvar_res.ok:
                                clinvar_json=clinvar_res.json()
                                variation_id=clinvar_json.get("primary_snapshot_data",{}).get("variation_id","")
                                if variation_id:
                                    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"
                                    clinvar_results.append({
                                        "id": var_id,
                                        "clinvar_link": clinvar_url,
                                    
                           
                                    })
                                    print(clinvar_results)

   
    return render_template("result.html",enter=gene_id,sequence=sequence,chromosome=chromosome,start=start,end=end,strand=strand, vep_results=vep_results,
                                clinvar_results=clinvar_results)
        ##here even though i didnt call element from result and still passed value as to pass value in html u just need to use render_template and make instance or element in the html

@app.route('/gene2', methods=['POST'])
def get_gene2():
    sequence2 = request.form.get('enter2','').strip()
    if not sequence2:
        return render_template("error.html",message="please enter a DNA sequence")
    
    if sequence2.startswith(">"):
            lines = sequence2.splitlines()
            if len(lines) == 1:
            # Handle pasted '>header sequence'
                parts = lines[0].split(None, 1)
                if len(parts) == 2:
                    fasta_seq = f"{parts[0]}\n{parts[1]}"
                else:
                    return render_template("error.html", message="No sequence data found after FASTA header.")
            else:
                fasta_seq = sequence2
    else:
        fasta_seq = f">Sequence\n{sequence2}"
        

    print(f"fasta submitted: {repr(fasta_seq)}")

    # Submit BLAST
    blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    submit_params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "nt",
        "QUERY": fasta_seq
    }

    submit_res = requests.post(blast_url, data=submit_params)
    print(f"fasta submitted: {repr(fasta_seq)}")
    clean_fasta="".join(line.strip() for line in fasta_seq if not line.startswith(">"))

    m=re.search(r"RID = ([\w-]+)",submit_res.text)
    if not m:
        print(submit_res.text)
        return "BLAST submission failed"
    rid=m.group(1)

    # Poll for results
    attempts=0
    max_attempts=60
    while attempts<max_attempts:
        check_params={
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo"
        }
        status_res = requests.get(blast_url, params=check_params)
        if "Status=READY" in status_res.text:
            break
        attempts+=1
        time.sleep(5)
        if attempts==max_attempts:
            return "Blast timedout"

    # Retrieve final result
    result_params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML"
    }
    xml_res= requests.get(blast_url, params=result_params)
    blast_record=NCBIXML.read(io.StringIO(xml_res.text))
     
    query_len=len(clean_fasta)
    top_hit=None
    top_score=float('-inf')
    mutate=None
    coverage=None
    # vep_results2=[]
    

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.score > top_score:
                top_score=hsp.score
                top_hit={
                     "title": alignment.title,
                    "score": hsp.score,
                    "evalue": hsp.expect,
                    "Identity":hsp.identities,
                    "query": hsp.query,
                    "match": hsp.match,
                    "subject": hsp.sbjct,
                    "start": hsp.query_start,
                    "end": hsp.query_end,
                    "start_sub":hsp.sbjct_start,
                    "end_sub":hsp.sbjct_end
                    
                }
                
                
    if top_hit:
        mutate=detect_mutation(hsp.query,hsp.sbjct)
        start = top_hit['start']
        end = top_hit['end']
        coverage = ( end - start + 1)/ query_len*100

        query_seq = top_hit['query']
        subject_seq = top_hit['subject']
        match_line = top_hit['match']


        q_start = top_hit['start']
        s_start = top_hit['start_sub']

        q_pos = q_start
        s_pos = s_start

        lines = []
        line_len = 60

        for i in range(0, len(query_seq), line_len):
            q_chunk = query_seq[i:i + line_len]
            m_chunk = match_line[i:i + line_len]
            s_chunk = subject_seq[i:i + line_len]

            q_end = q_pos + len(q_chunk) - 1
            s_end = s_pos + len(s_chunk) - 1

            lines.append(f"Query   {q_pos:<5} {q_chunk} {q_end}")
            lines.append(f"        {'':<5} {m_chunk}")
            lines.append(f"Subject {s_pos:<5} {s_chunk} {s_end}\n")

            q_pos += len(q_chunk)
            s_pos += len(s_chunk)

# Join formatted alignment
        formatted_alignment = "\n".join(lines)

        
        top_result=(
            f"Title: {top_hit['title']}\n\n"
            f"Score: {top_hit['score']}\n"
            f"E-value: {top_hit['evalue']}\n"
            f"Identity:{top_hit['Identity']}\n"
            f"Query length:{query_len}\n"
            f"Query Coverage:{coverage:.2f}%\n\n"
             f"{formatted_alignment}"
             
           
        )
    else:
        return "No top hits found"

    
    return render_template("ukresult.html", enter2=sequence2, blast_output=top_result,Mutagen=mutate)

if __name__ == '__main__':
    app.run(debug=True)