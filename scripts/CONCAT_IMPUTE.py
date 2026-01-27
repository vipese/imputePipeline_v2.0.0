"""
scripts to parse impute2 output and concat the imputed genotypes and concat the files
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""
import re
import sys
import argparse
import gzip
import datetime
import os
import glob

dttime=datetime.datetime.now().strftime ("%Y%m%d")
parser = argparse.ArgumentParser(description='A class method to parse impute2 output and concat the imputed genotypes and concat the files')
parser.add_argument('-F', required=True, help='Suffix of the plink binaries used in the head command for impute2')

### parse the arguments
args=parser.parse_args()
file_prefix = args.F

### make a class object to clean up the folder ####
class impute_clean(object):
    """A class of impute object:

    Attributes:
        name : A string identifier for the object class
        file_prefix : A file prefix to be used to parse the directory

    """
    instance_count =0

    def __init__(self, name, file_prefix):
        """intitate the LDmap object class."""
        self.name = name
        self.file_prefix = file_prefix
        self.instance_count += 1

    def get_attr(self):
        print(' NAME :- {} \n SOURCE FILE PREFIX:- {} \n INSTANCE COUNT :- {}'.format(self.name, self.file_prefix,  self.instance_count))

    def parse_dir(self):
        file_ = self.file_prefix
        getcwd =os.getcwd()
        imp_accfiles = {'warnings':'', 'summary':'', 'info':'', 'info_by_sample':'', 'diplotype_ordering':'', 'impute':''}
        for filename in os.listdir(getcwd):
            #print filename
            if file_ in filename:
                assert file_ in filename
                basename, ext = os.path.splitext(filename)
                #print(basename, ext)
                if ext.replace('.','').isdigit():
                    make_key = 'impute'
                    #print 'yeay'
                else:
                    make_key = '_'.join(ext.split('_')[1:])
                if make_key in imp_accfiles:
                    get_exist = imp_accfiles.get(make_key)
                    imp_accfiles[make_key]= get_exist+'\n'+filename
                else:
                    pass
        return imp_accfiles
        
    @staticmethod
    def write_out(filestring, outfile):
        if os.path.exists(filestring):
            with open(filestring) as infile:
                line_buf =0
                line_n =0
                for line in infile:
                    line_n += 1
                    if line_n == 10000:
                        line_buf += 10000
                        line_n =0
                        print('PROCESSED {} LINES '.format(line_buf))
                    outfile.write(line)
                #outfile.close()

    #@staticmethod
    def concat_files(self):
        file_dic=self.parse_dir()
        if len(file_dic) == 6:
            assert len(file_dic) == 6
            for k, v in file_dic.items():
                make_outfile = self.file_prefix+'.'+k+'.gz'
                print('MADE A GZIP FILE CALLED {} '.format(make_outfile))
                if os.path.exists(make_outfile):
                    outfile = gzip.open(make_outfile, 'a')
                    print('FILE EXISITS APPENDING TO ARCHIVE {} '.format(outfile))
                else:
                    outfile = gzip.open(make_outfile, 'wb')
                parse_files = v.split('\n')
                for filestring in parse_files:
                    self.write_out(outfile=outfile, filestring=filestring)




if __name__ == '__main__':
    impute_object = impute_clean(name='myimpute', file_prefix=file_prefix)
    impute_object.get_attr()
    #impute_object.parse_dir()
    impute_object.concat_files()
#file_='CHR6_Plates_77_to_109_PMRA_shapeit'#124
 