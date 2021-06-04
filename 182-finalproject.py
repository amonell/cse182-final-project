#!/usr/bin/env python
# coding: utf-8

# In[13]:


import gzip
with gzip.open('GCF_009858895.2_ASM985889v3_genomic.gff.gz', 'rb') as f:
    file_content = str(f.read())


# In[14]:


splitted = file_content.split("\\n")
splitted = [i for i in splitted if i[0] == "N" and i[1] == "C"]


# In[15]:


first = splitted[0]
first
first = first.split("\\t")
#for i in first:
#    print(i)


# In[16]:


proteins = [j for j in splitted if "protein_id" in j]
headers = [j.split("protein_id")[0] for j in proteins]
proteins = [j.split("protein_id")[1][1:] for j in proteins]


# In[17]:


x = None
with open("AccessionNumbers.txt") as file:
    x = file.readlines()
    x = [i for i in x if ">" not in i and i != "\n"]
    for j in range(len(x)):
        if "\n" in x[j]:
            x[j] = x[j][:len(x[j])-1]
number = 1 + ((2-1)*400)%10000
accesions = x[number:number+150]
#accesions


# In[18]:


query = None
with open("nucseq.fasta.txt") as file:
    query = file.readlines()
    for j in range(len(query)):
        if "\n" in query[j]:
            query[j] = query[j][:len(query[j])-1]

query_dict = {}
indexes = []
count = 0
for i in query:
    #print(i)
    if len(i) > 0:
        if i[0] == ">":
            indexes.append(count)

    count += 1

def Tostring(sublist):
    stri = ""
    for i in sublist:
        stri += i
    return stri

for i in range(len(indexes)):
    #print(i)
    sublist = None
    if i != len(indexes) - 1:  
        sublist = query[indexes[i]+1:indexes[i+1]]
    else:
        sublist = query[indexes[i]+1:]
        
    query_dict[query[indexes[i]]] = Tostring(sublist)
    

#print(query_dict['>MT502956.1 Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/THA/CO63-896/2020 ORF1ab polyprotein (ORF1ab) and ORF1a polyprotein (ORF1ab) genes, partial cds'])


# In[19]:



def Tostring(sublist):
    stri = ""
    for i in sublist:
        stri += i
    return stri

protein_database = None
with open("proteinsequences.fasta.txt") as file:
    protein_database = file.readlines()
    for j in range(len(protein_database)):
        if "\n" in protein_database[j]:
            protein_database[j] = protein_database[j][:len(protein_database[j])-1]

protein_dict = {}
indexes = []
count = 0
for i in protein_database:
    #print(i)
    if len(i) > 0:
        if i[0] == ">":
            indexes.append(count)

    count += 1

def Tostring(sublist):
    stri = ""
    for i in sublist:
        stri += i
    return stri

for i in range(len(indexes)):
    #print(i)
    sublist = None
    if i != len(indexes) - 1:  
        sublist = protein_database[indexes[i]+1:indexes[i+1]]
    else:
        sublist = protein_database[indexes[i]+1:]
        
    protein_dict[protein_database[indexes[i]]] = Tostring(sublist)


# In[20]:


newest_dic = {}
for i in protein_dict:
    newest_dic[i] = len(protein_dict.get(i))


# In[21]:


reference_genome = None
with open("covidgenome.fasta.txt") as file:
    reference_genome = file.readlines()
    for j in range(len(reference_genome)):
        if "\n" in reference_genome[j]:
            reference_genome[j] = reference_genome[j][:len(reference_genome[j])-1]
reference_genome
string = Tostring(reference_genome[1:])
reference_genome = {reference_genome[0]: string}
#reference_genome


# In[ ]:


#blast output has each individual query and if it has any matches with the protein database.
#We want to know in each individual protein where there is a difference between the actual protein versus what our query sequence has.
#


# In[22]:


blast = None
with open("outputfileblast2.txt") as file:
    blast = file.readlines()
    for j in range(len(blast)):
        if "\n" in blast[j]:
            blast[j] = blast[j][:len(blast[j])-1]

query_indexes = []
count = 0
for i in blast:
    #print(i)
    if len(i) > 0:
        if "Query=" in i:
            query_indexes.append(count)

    count += 1

test_set = blast[query_indexes[1]:query_indexes[2]]

#test_set


# In[23]:


def find_mutations(test_set, protein_mutations):
    query_name = test_set[0].split()[1]

    query_length = None
    start_bitscores = None

    print("Query = " + query_name)
                       
    for i in range(len(test_set)):
        if "Length" in test_set[i]:
            query_length = test_set[i].split("=")[1]
            start_bitscores = i + 4
            break

    protein_table = {}
    for i in range(start_bitscores, len(test_set)):
        if test_set[i][0:2] == "YP":
            if float(test_set[i].split()[len(test_set[i].split())-2]) >= 50.0:
                protein_table[test_set[i].split()[0]] = test_set[i].split()[len(test_set[i].split())-2]
        else:
            break

    protein_indexes = []
    count = start_bitscores
    for i in test_set[start_bitscores:]:
        #print(i)
        if len(i) > 0:
            if i[0] == ">":
                protein_indexes.append((count, i))

        count += 1

    sections = {}
    for i in range(len(protein_indexes)):
        #Collect only query lines and the 2 lines after
        if i != len(protein_indexes) - 1:
            sections[protein_indexes[i][1].split()[0]] = test_set[int(protein_indexes[i][0]):int(protein_indexes[i+1][0])]
        else:
            sections[protein_indexes[i][1].split()[0]] = test_set[int(protein_indexes[i][0]):]

    for i in sections:
        framesections = -1
        frame_ct = 0
        for j in sections.get(i):
            if "Frame" in j:
                framesections = frame_ct
                break
            frame_ct += 1
        sections[i] = sections.get(i)[framesections + 2:]

    for i in sections:
        framesections = len(sections.get(i))
        frame_ct = 0
        for j in sections.get(i):
            if "Lambda" in j:
                framesections = frame_ct
                break
            frame_ct += 1
        sections[i] = sections.get(i)[:framesections-1]

    new_additions = {}
    for i in sections:
        frame_ind = []
        frame_ct = 0
        for j in sections.get(i):
            if "Frame" in j:
                frame_ind.append(frame_ct)
            frame_ct += 1
        indcount = len(frame_ind)
        for k in range(indcount):
            if k == 0:
                if k != indcount-1:
                    new_additions[i+str(k)] = sections.get(i)[frame_ind[k]+1:frame_ind[k+1]]
                    new_additions[i] = sections.get(i)[:frame_ind[k]]
                else:
                    new_additions[i+str(k)] = sections.get(i)[frame_ind[k]+1:]     
                    new_additions[i] = sections.get(i)[:frame_ind[k]]
            if k != indcount-1:
                new_additions[i+str(k)] = sections.get(i)[frame_ind[k]+1:frame_ind[k+1]]
            else:
                new_additions[i+str(k)] = sections.get(i)[frame_ind[k]+1:] 

    for i in new_additions:
        sections[i] = new_additions.get(i)


    for i in sections:
        jcount = 0
        for j in sections.get(i):
            if "Score" in j:
                sections[i] = sections.get(i)[:jcount]
            jcount += 1
            
    for i in sections:
        sections[i] = [k for k in sections.get(i) if len(k) > 0]
    for i in sections:
        new_setting = []
        for j in range(2, len(sections.get(i)), 3):
            #print(j, len(sections.get(i)))
            new_setting.append([sections.get(i)[j-2], sections.get(i)[j-1], sections.get(i)[j]])
        sections[i] = new_setting

    query_set = set([i for i in range(int(query_length))])
    
    for i in sections:
        if i not in protein_mutations:
            protein_mutations[i] = set({})
            
    for i in sections:
        counterj = 0
        for j in sections.get(i):
            if counterj == 0:
                if j[0].split()[2].lower() != j[2].split()[2].lower():

                    starter = int(j[2].split()[1])
                    og_starter = starter
                    query_start = int(j[0].split()[1])
                    query_end = int(j[0].split()[3])
                    nuc_splitter = j[0].split()[2].lower()
                    splitter = j[2].split()[2].lower()
                    for let in range(len(splitter)):
                        if splitter[let].lower() != nuc_splitter[let].lower():
                            if starter != og_starter and query_start in query_set:
                                print(i + ":" + splitter[let].lower() + str(starter) + nuc_splitter[let].lower())
                                protein_mutations[i].add(i + ":" + splitter[let].lower() + str(starter) + nuc_splitter[let].lower())
                                #query_set.remove(query_start)
                            elif starter == og_starter and query_start in query_set:
                                print(i + ":" + "indel" + str(starter) + "indel")
                                protein_mutations[i].add(i + ":" + "indel" + str(starter) + "indel")
                                og_starter += 1
                            elif starter == og_starter and query_start not in query_set:
                                og_starter += 1                                
                                
                        starter += 1
                        query_start += 1

            elif j[0].split()[2].lower() != j[2].split()[2].lower():
                #print(i)
                #print(j[0].split()[2])
                #print(j[2].split()[1])
                #print(j[2].split()[2])
                #print(j[2].split()[3])
                starter = int(j[2].split()[1])
                query_start = int(j[0].split()[1])
                query_end = int(j[0].split()[3])
                nuc_splitter = j[0].split()[2].lower()
                splitter = j[2].split()[2].lower()
                for let in range(len(splitter)):
                    if splitter[let].lower() != nuc_splitter[let].lower() and query_start in query_set:
                        print(i + ":" + splitter[let].lower() + str(starter) + nuc_splitter[let].lower())
                        protein_mutations[i].add(i + ":" + splitter[let].lower() + str(starter) + nuc_splitter[let].lower())
                        #query_set.remove(query_start)
                    starter += 1
                    query_start += 1
            counterj += 1


# In[24]:


blast = None
with open("outputfileblast2.txt") as file:
    blast = file.readlines()
    for j in range(len(blast)):
        if "\n" in blast[j]:
            blast[j] = blast[j][:len(blast[j])-1]

query_indexes = []
count = 0
for i in blast:
    #print(i)
    if len(i) > 0:
        if "Query=" in i:
            query_indexes.append(count)

    count += 1

protein_mutations = {}

for i in range(1, int(len(query_indexes)/2+1)):
    test_set = blast[query_indexes[i-1]:query_indexes[i]]
    find_mutations(test_set, protein_mutations)

print("Mutations per Protein")
for i in protein_mutations:
    print(i)
    print("-----------------------------------")
    for j in protein_mutations.get(i):
        print(j)
    print("-------------------------------------")

