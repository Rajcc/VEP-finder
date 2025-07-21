def detect_mutation(query,subject):
    mutations=[]

    complement = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C'
    }
    purines=['A','G']
    pyrimidines=['C','T']

    for i in range(len(query)):
        q_char=query[i]
        s_char=subject[i]
        kind=None

        if q_char==s_char:
            continue
        elif q_char=='-':
            mutations.append({'Type':'Deletion','Position':i+1,'Query Base':'-','Subject Base':s_char,'Variation kind':'Indel'})
        elif s_char=='-':
            mutations.append({'Type':'Insertion','Position':i+1,'Query Base':q_char,'Subject Base':'-','Variation kind':'Indel'})
        else:
            if complement.get(q_char)!=s_char:
                if (q_char in purines and s_char in purines) or (q_char in pyrimidines and s_char in pyrimidines):
                    kind='Transition(SNP)'
                else:
                    kind='Transversion(SNP)'

                mutations.append({'Type':'MisMatch','Position':i+1,'Query Base':q_char,'Subject Base':s_char,'Variation kind':kind})
    
    return mutations