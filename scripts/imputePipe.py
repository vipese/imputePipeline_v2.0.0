import re
from subprocess import Popen, PIPE
import sys
import os
import argparse


def checkPath(filePre, fileSuf, chr_):
    '''Given a prefix and file extension and a chr number checks if a file exists in current directory
    '''
    pathFile = os.path.join(os.getcwd(), '{}_CHR{}.{}'.format(filePre, chr_, fileSuf))
    print(pathFile)
    if os.path.exists(pathFile):
        print('FOUND MATCHING FILE {}'.format(pathFile))
        return True
    else:
        return False


def plinkSplitCall(filePre):
    '''Given a base bed file, this function will submit a job to split files by chromosomes and return a job ID
    '''
    if os.path.exists(os.path.join(os.getcwd(), 'scripts/PLINK_SPLIT_SLURM.sh')):
        plinKSplit='./scripts/PLINK_SPLIT_SLURM.sh {}'.format(filePre)
        plinkCall=Popen(plinKSplit, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr= plinkCall.communicate()
        if not stderr:
            job1 = re.findall(r'\d+', stdout.decode())[0]
            if not job1:
                exit('Job Submission Unsuccessful for Script {}'.format(plinKSplit))
            else:
                return job1
        else:
            exit('Job Submission Unsuccessful for Script {}'.format(plinKSplit))
    else:
        raise FileNotFoundError('Check Path scripts/PLINK_SPLIT_SLURM.sh')



def phasingCall(filePre, *args):
    '''function to take input as the base plink file prefix and an optional job dependency for slurm
    '''
    if args:
        dependency=args[0]
        if dependency:
            print('Job Dependency {}'.format(dependency))
            shapeIt='./scripts/SHAPEIT_ARRAY_TASK_SLURM.sh {} {}'.format(filePre, dependency)
    else:
        shapeIt='./scripts/SHAPEIT_ARRAY_TASK_SLURM.sh {}'.format(filePre)
    print(shapeIt)
    if os.path.exists(os.path.join(os.getcwd(), 'scripts/SHAPEIT_ARRAY_TASK_SLURM.sh')):
        shapeitCall=Popen(shapeIt, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr= shapeitCall.communicate()
        if not stderr:
            job2 = re.findall(r'\d+', stdout.decode())[0]
            if not job2:
                exit('Job Submission Unsuccessful for Script {}'.format(shapeitCall))
            else:
                return job2
        else:
            exit('Job Submission Unsuccessful for Script {}'.format(shapeitCall))
    else:
        raise FileNotFoundError('Check Path scripts/SHAPEIT_ARRAY_TASK_SLURM.sh')


def imputeCall(filePre, ref, *args):
    '''makes the imputation calls with job dependency if applicable to the shapeit call
    '''
    chrSizes = [('1', '249250621'), 
                ('2', '243199373'), 
                ('3', '198022430'), 
                ('4', '191154276'), 
                ('5', '180915260'), 
                ('6', '171115067'), 
                ('7', '159138663'), 
                ('8', '146364022'), 
                ('9', '141213431'), 
                ('10', '135534747'), 
                ('11', '135006516'), 
                ('12', '133851895'), 
                ('13', '115169878'), 
                ('14', '107349540'), 
                ('15', '102531392'), 
                ('16', '90354753'), 
                ('17', '81195210'), 
                ('18', '78077248'), 
                ('20', '63025520'), 
                ('19', '59128983'),
                ('21', '48129895'),
                ('22', '51304566')]
    if ref == "3":
        imputeScript = 'IMPUTE_LOOP_SLURM_FIXED.sh'
    else:
        imputeScript = 'IMPUTE_LOOP_SLURM_Phase1.sh'

    if args:
        dependency=args[0]
        if dependency:
            print('Job Dependency {}'.format(dependency))
            for chrTup in chrSizes:
                chr_=chrTup[0]
                chrLength=chrTup[1]
                if len(chrLength) < 9:
                    chrSize = str(int(chrLength[:2])+1)
                else:
                    chrSize = str(int(chrLength[:3])+1)
                impute='./scripts/{} {} {} {} {}'.format(imputeScript, chr_, chrSize, filePre, dependency)
                Popen(impute, shell=True, stdout=PIPE, stderr=PIPE)
                print(impute)
    else:
        for chrTup in chrSizes:
            chr_=chrTup[0]
            chrLength=chrTup[1]
            if len(chrLength) < 9:
                chrSize = str(int(chrLength[:2])+1)
            else:
                chrSize = str(int(chrLength[:3])+1)
            impute='./scripts/{} {} {} {}'.format(imputeScript, chr_, chrSize, filePre)
            Popen(impute, shell=True, stdout=PIPE, stderr=PIPE)
            print(impute)
            
def main():
    parser = argparse.ArgumentParser(description='Imputation Pipeline Main')
    parser.add_argument('-F', help='File Prefix for the base BED file unspit by chromosomes', required=True)
    parser.add_argument('-Ref', help='1000 genomes reference the input should be either 1 or 3', required=True)
    args=parser.parse_args()
    filePre=args.F
    ref=args.Ref
    plinkFile=checkPath(filePre, fileSuf='bed', chr_='2')
    hapFile=checkPath(filePre, fileSuf='haps', chr_='2')
    if plinkFile: ## if plinkfile exists check for hap file
        if hapFile: ## if hap file exists go to imputation
            imputeCall(filePre, ref)
        else:
            job2=phasingCall(filePre)
            if job2:
                imputeCall(filePre, ref, job2)    
    else:
        job1 = plinkSplitCall(filePre)
        if job1:
            job2 = phasingCall(filePre, job1)
            if job2:
                imputeCall(filePre, ref, job2)

    

if __name__ == "__main__":main()