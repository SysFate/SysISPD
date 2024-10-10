import argparse
from xml.dom.minidom import parseString
from xmlrpc.client import Boolean
import pandas as pd
import sys
import os
import shutil
import pysam

from os import path,listdir
from datetime import datetime
from subprocess import Popen, PIPE

ID = None
tmp = "/tmp"
tsvFolder = "./"
cmdDecompress = ["pigz", "-dc"]
core = 1
featuresToGenerate = ["exon","transcript"]

def quantification(samDir, feature, gtfFile, dictConvertPositionMatrix, cmdQuantification):
    sys.stdout.write(f'### - Quantification - {feature}\n')
    cmd_count = []
    global tsvFolder
    pathTSVFile = path.join(tsvFolder,f'matrix_count_{feature}.tsv')
    #cmdQuantification = "featureCounts -T __core__ -F GTF -p --countReadPairs -t __feature__ -s 0 -a __GTF__ -o __out__"
    for partCMD in cmdQuantification.split(' '):
        if partCMD == '__GTF__':
            cmd_count.append(gtfFile)
        elif partCMD == '__core__':
            cmd_count.append(str(core))
        elif partCMD == '__feature__':
            cmd_count.append(feature)
        elif partCMD == '__out__':
            cmd_count.append(pathTSVFile)
        else:
            cmd_count.append(partCMD)
    dictConvertsamFileBarcode = {}
    for filename in listdir(samDir):
        file = path.join(samDir, filename)
        if path.isfile(file) and filename.endswith('.sam') and path.getsize(file)>0:
            barcode=filename.strip('.sam')
            cmd_count.append(file)
            dictConvertsamFileBarcode[file] = barcode
   # for barcode in listBarcode:
        #samFile = path.join(samDir, f'{barcode}.sam')
        #cmd_count.append(samFile)
        #dictConvertsamFileBarcode[samFile] = barcode

    #for element in cmd_count:
        #print(element)

    quantiRun = Popen(cmd_count, stdout=PIPE, stderr=PIPE)
    stdout, stderr = quantiRun.communicate()
    with open(path.join(tmp,f"quantification_{feature}.log"), 'w') as log, open(path.join(tmp,f"quantification_{feature}.err"), 'w') as error:
        log.write(stdout.decode('utf-8'))
        error.write(stderr.decode('utf-8'))
    if quantiRun.returncode != 0:
        raise Exception(f'Error with counting {feature}')
    #abs_file_path = '/home/mduvina/Documents/Cnrs/python/epigenomic/sam_files/'
    #pathTSVFile = path.join(abs_file_path,'output.tsv')
    countMatrix = pd.read_csv(pathTSVFile, sep = '\t', skiprows=1, index_col=0) # <= skiprows is for the first row of feature coun
    summaryCount = pd.read_csv(pathTSVFile+".summary", sep = '\t', index_col=0, header=0)
    #### Clean special feature count ####
    
    for col in ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']:
        try:
            countMatrix.drop(col, axis=1, inplace=True)
        except:
            pass
    countMatrix = countMatrix[~(countMatrix == 0).all(axis=1)]
    countMatrix.rename(dictConvertsamFileBarcode, axis=1, inplace=True)
    countMatrix.rename(dictConvertPositionMatrix, axis=1, inplace=True)
    totalCount = int(countMatrix.to_numpy().sum())
    countMatrix.to_csv(pathTSVFile, sep='\t')
    summaryCount.rename(dictConvertsamFileBarcode, axis=1, inplace=True)
    summaryCount.rename(dictConvertPositionMatrix, axis=1, inplace=True)
    summaryCount.to_csv(pathTSVFile+".summary", sep='\t')
    
    print(totalCount)

def editQuantifierCommand(quantifier, ignoreDup):
    if ignoreDup:
        quantifier = quantifier+" --ignoreDup"
    print(quantifier)
    return quantifier

    
def make_dirs(d):
    try:
        os.makedirs(d)
        return True
    except:
        return False
    
def test_file_can_open(file,errorMessage):
    try:
        with open(file):
            pass
    except IOError:
        raise Exception(errorMessage)

def remove_dir(dir_path,printError=True):
    try:
        shutil.rmtree(dir_path)
    except OSError as e:
        if printError:
            print("Error: %s : %s" % (dir_path, e.strerror))
        
def move_dir(Dir,Folder):
    try:
        shutil.move(Dir,Folder)
        return True,None
    except Exception as e:
        return False,e

def extract_barcode_present(matrix_position):
    dictConvertPositionMatrix = {}
    df = pd.read_csv(matrix_position,sep="\t",index_col=0)
    for index, row in df.iterrows():
        for nameCol in list(df):
            if isinstance(row[nameCol],str):
                results = row[nameCol][1:].split('y')
                dictConvertPositionMatrix[results[0]+'x'+results[1]] = '{}x{}'.format(index,nameCol)
    return dictConvertPositionMatrix
    
def modGTFfile(gtfFile, modgtf):
    if modgtf == 0:
        return gtfFile
    print("##### Changing GTF #####")
    name = path.basename(gtfFile).strip(".gtf")
    oldGtf = open(gtfFile)
    newGtf = open(path.join(tmp,name+f"+-{modgtf}kb.gtf"), "a")
    for line in oldGtf:
        line=line.split("\t")
        if line[6] == '+':
            line[4]=str(int(line[3])+modgtf*1000)
            line[3]=str(max([int(line[3])-modgtf*1000,1]))
        else:
            line[3]=str(max([int(line[4])-modgtf*1000,1]))
            line[4]=str(int(line[4])+modgtf*1000)
        line="\t".join(line)
        newGtf.write(line)
    return path.join(tmp,name+f"+-{modgtf}kb.gtf")


def mergeSam(samDir):
    list_of_files = []
    for filename in listdir(samDir):
        file = path.join(samDir, filename)
        if path.isfile(file) and filename.endswith('.sam') and path.getsize(file)>0:
            list_of_files.append(file)
    output=path.join(tmp,"bulk_sam.sam")
    pysam.merge("-f","-o",output,*list_of_files)
    

def main():
    global ID
    global outFolder
    global tsvFolder
    global featuresToGenerate
    global tmp
    global cmdDecompress
    global core
    global cmdQuantification
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n","--name", metavar='name', type=str, default=datetime.now().strftime("%Y-%m-%d_%H:%M:%S"), help="You can give a name for your analyse. By default the software use the begin date.")
    parser.add_argument("--tmp", metavar='tmp', type=str, default="./tmp", help="Folder where send temporary files")
    parser.add_argument("-o","--out", metavar='outdir', type=str, default=".", help="Folder where send files")
    parser.add_argument("--decompress", metavar='decompress', type=str, default="pigz -dc", help="Command use to decompress fastq files")
    #Importants:
    parser.add_argument("-q","--quantifier", metavar='quantifier', type=str, default="featureCounts -T __core__ -F GTF -t __feature__ -s 0 -a __GTF__ -o __out__ -O", help="This is the read quantifier command. You can write what you want if the __GTF__, __feature__ is provide and the output is configured to send the results in __out__.")
    parser.add_argument("-c","--core", metavar='core', type=int, default=1, help="Number of core the software can use")
    parser.add_argument("--feature", metavar='feature', type=str, default="5UTR,transcript", help="Feature use during the read quantifier in your file separate by a comma")
    parser.add_argument("--gtf", metavar='gtf', type=str, required=True, help="GTF file to use with feature count")
    parser.add_argument("--pos", metavar='position', type=str, required=True, help="TSV file who contain coordinates code with their position in the grid")
    #New ones:
    parser.add_argument("--modgtf", metavar='modgtf', type=int, default=0, help ="Number of kb to ajust GTF")
    parser.add_argument("--samdir", metavar="samdir", type=str, required=True, help ="Directory of SAM files")
    parser.add_argument("--merge", action = 'store_true', help ="Merge sam files")
    parser.add_argument("--ignoreDup", action = 'store_true', help = "Ignore duplicate")
    args = parser.parse_args()
    test_file_can_open(args.gtf, "GTF file can't be open")
    ID = args.name
    matrix_position = args.pos
    outFolder = path.join(args.out,args.name)
    if path.exists(outFolder):
        i = 0
        while path.exists(outFolder+"."+str(i)):
            i = i+1
        outFolder = outFolder+"."+str(i)
    make_dirs(outFolder)
    sys.stdout.write(f'### - {args.name}\n')
    tmp = path.join(args.tmp, f'{args.name}_demultiplex_tmp')
    if path.exists(tmp):
        i = 0
        while path.exists(tmp+"."+str(i)):
            i = i+1
        tmp = tmp+"."+str(i)
    make_dirs(tmp)
    cmdDecompress = args.decompress.split(" ")
    core = args.core
    samDir = args.samdir
    gtfFile = args.gtf
    modgtf = args.modgtf
    gtfFile = modGTFfile(gtfFile, modgtf)
    featuresToGenerate = args.feature.split(",")
    tsvFolder = outFolder
    merge = args.merge
    #python3 test.py --samdir ../../sam_files/ --gtf ../../sam_files/mm9.gtf --pos ../../sam_files/Interstitial_Matrix_32x64_V2.tsv 
    #--modgtf 0 --quantifier "featureCounts -T __core__ -F GTF -p --countReadPairs -t __feature__ -s 0 -a __GTF__ -o __out__"
    #samDir='/home/mduvina/Documents/Cnrs/python/epigenomic/sam_files/' #à faire (juste demander l'argument)
    dictConvertPositionMatrix = extract_barcode_present(matrix_position)
    if merge:
        mergeSam(samDir)
    else:
        cmdQuantification = editQuantifierCommand(args.quantifier, args.ignoreDup)
        for feature in featuresToGenerate:
            quantification(samDir, feature, gtfFile, dictConvertPositionMatrix, cmdQuantification)   
    #remove_dir(tmp)

if __name__ == '__main__':
    main()