from os import listdir, rename
from os.path import isfile, join, isdir
#import multiprocessing #추후 멀티프로세싱 도입해야할듯
import numpy as np

gene_path =  "./gene/geneFamily"
familyNames = []
Kmer = 20

def makeFamilyFragmentDictionaries():#각 Family별 gene Dictionary 생성
    familyNames.extend( [f for f in listdir(gene_path) if isdir(join(gene_path, f))] )#전역변수 gene_path 폴더안에 있는 폴더 검색
    paths = [join(gene_path, f) for f in familyNames]
    sequenceFragmentForEachFamily = [] #family별 단일 sequence의 fragment set
    dictionariesForEachFamily = [] #family별 전체 fragment dictionary

    for path in paths:#paths에는 family directory folder 경로
        sequences = makeFamilySequences(path)
        dictionary, fragmentForEachSequence = makeFragmentFrequencyDict(sequences)
        sequenceFragmentForEachFamily.append(fragmentForEachSequence)
        dictionary = {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1], reverse=True)}#fragment 등장 빈도를 기준으로 sorting
        dictionariesForEachFamily.append(dictionary)
    #print(sequenceFragmentForEachFamily)
    return dictionariesForEachFamily, sequenceFragmentForEachFamily

#makeFamilyDictionaries() 에서 사용
def makeFamilySequences(family_dir_loc):#각 Family별 sequence 리스트 생성
    filenames = [f for f in listdir(family_dir_loc) if isfile(join(family_dir_loc, f))]#각 family 폴더안에 있는 Fasta파일 검색
    sequences = []
    for filename in filenames:
        path = join(family_dir_loc, filename)
        with open(path, "r") as f:
            head, tail = f.read().split('\n', 1)#Fasta format 문서 첫줄제거
            sequence = tail.replace("\n", "")#한줄로
            sequences.append(sequence)
    return sequences

#makeFamilyDictionaries() 에서 사용
def makeFragmentFrequencyDict(sequences):#sequence를 kmer사이즈로 잘라서 사전에 추가, 사전에는 fragment의 발생 빈도 기록
    fragmentForEachSequence = []#각각의 gene에 속해있는 fragment
    dictionary = dict()#사전 생성
    for sequence in sequences:
        length = len(sequence)
        fragments = set()#중복없는 set으로 fragment
        for i in range(0, length - Kmer + 1):
            frag = sequence[i: i + Kmer]#k개씩 자름
            fragments.add(frag)

        for fragment in fragments:
            if fragment in dictionary:
                dictionary[fragment] += 1#이미 있으면 +1
            else:
                dictionary[fragment] = 1#없으면 1

        fragmentForEachSequence.append(fragments)
    #print(fragmentForEachSequence)
    return dictionary, fragmentForEachSequence  #dictionary - key: fragment, value: fragment 발생 빈도,  fragmentForEachSequence - [{sequence별 fragment}, {} ...]

def makeOverlappingSeqeunceResult(dictionariesForEachFamily):#패밀리별 fragment를 입력받아 각 fragment가 각각의 패밀리에 속해 있는지 검사
    result = []
    #모든 family의 fragment가 속한 집합(fullSet) 생성
    fullSet = set()
    for family_dict in dictionariesForEachFamily:
        key = family_dict.keys()
        fullSet.update(key)

    #fullSet에 속한 fragment들을 각각의 Family에 속해있는지 검사
    for setItem in fullSet:
        line = setItem
        sum = 0

        for family_dict in dictionariesForEachFamily:
            if setItem in family_dict:
                freq = family_dict[setItem]
                line += "," + str(freq)
                sum += freq
            else:
                line += "," + str(0)
        line += "," + str(sum)
        result.append([line, sum])
    result = sorted(result, key=lambda item: item[1], reverse=True)
    return result

def makeFragmentFrequencyDictToCSV(fragmentFrequencies):
    for i in range(len(fragmentFrequencies)):
        with open("./fragment analysis/" + familyNames[i] + "_FragmentFrequencyResult"+"_Kmer="+str(Kmer) +".csv", "w") as f:
            f.write("fragment,frequency\n")
            for key,value in fragmentFrequencies[i].items():
                f.write(key + "," + str(value) +"\n")

def makeOverlappingResultToCSV(overlappingResult):#familyNames: 이름을 나열한 리스트, fragmentFrequency: makeOverlappingSeqeunceResult의 리턴값
    with open("./fragment analysis/OverlappingResult"+"_kmer="+str(Kmer)+".csv", "w") as f:
        f.write("fragment")
        for familyName in familyNames:
            f.write("," + familyName)
        f.write(",sum\n")
        for element in overlappingResult:
            f.write(element[0] + "\n")

def convertFragmentToNumber(fragment):#Fragment를 숫자로 변형 A:1, G:2, T:3, C:4
    result = []
    for char in fragment:
        result.append(AGTCSwitcher(char))
    return result

def AGTCSwitcher(char):
    switcher={'A':1, 'G':2, 'T':3, 'C':4,
              "R":5, "Y":6, "K":7, "M":8}
    return switcher[char]

#clustering을 위한 함수
#전처리
# dictForEachFam:makeFamilyFragmentDictionaries() 결과 , frequency: 추출할 유전자 빈도수, 해당 빈도 이상의 fragment를 추출하기 위한 값
def preprocessingForClustering(dictForEachFam):
    x = []
    y = []
    for i in range(len(dictForEachFam)):
        max_frequency = max(dictForEachFam[i].values())
        frequency = max_frequency/2
        print("max_frequency = " + str(max_frequency) + ", frequency: " + str(frequency))

        for key, value in dictForEachFam[i].items():
            #if value >= 2:
            if value >= frequency and value > 1:
                x.append( convertFragmentToNumber(key) )
                y.append(i)
    x = np.array(x)
    y = np.array(y)
    #print(x)
    #print(y)
    return x,y

# fragment 발생 빈도 기록 dictionary 생성
print("Start process!")
print("makeFamilyFragmentDictionaries")
dictionariesForEachFamily, sequenceFragmentForEachFamily = makeFamilyFragmentDictionaries()
print("makeFragmentFrequencyDictToCSV")
makeFragmentFrequencyDictToCSV(dictionariesForEachFamily)

# 패밀리별 fragment를 입력받아 각 fragment가 각각의 패밀리에 속해 있는지 검사 후 결과 출력
print("makeOverlappingSeqeunceResult")
result = makeOverlappingSeqeunceResult(dictionariesForEachFamily)
print("makeOverlappingResultToCSV")
makeOverlappingResultToCSV(result)
print("Finished!")
# clustering을 위한 전처리
#x, y = preprocessingForClustering(dictionariesForEachFamily)

