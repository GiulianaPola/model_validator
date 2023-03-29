#!/usr/bin/env python3

version='1.0.8B'

import os
from datetime import datetime
import argparse
import sys
import subprocess
from Bio import Entrez
import time

param=dict()
validated=[]
invalidated=[]
error=[]
taxonomy=dict()
start_time = datetime.now()
call=os.path.abspath(os.getcwd())
finalresults=dict()

print('model_validator v{} - a tool to validate profile HMMs and adjust their cutoff score.\n'.format(version))

parser = argparse.ArgumentParser(prog = 'model_validator',add_help=False)
parser.add_argument('-i','-input_file')
parser.add_argument('-out','-output', '-o',default='hmm_validated')
parser.add_argument('-conf')
parser.add_argument('-db_type','-dbtype',type = str.lower,choices=['long', 'l', 'short','s'])
parser.add_argument('-model_type','-modeltype',type = str.lower,choices=['long', 'l', 'short','s'])
parser.add_argument('-cell_org_db')
parser.add_argument('-cell_org_db_frag')
parser.add_argument('-vir_db')
parser.add_argument('-pt',default=80)
parser.add_argument('-pd',default=80)
parser.add_argument('-long_hmms_x_short_reads_cell','-lh_sr')
parser.add_argument('-short_hmms_x_short_reads_cell','-sh_sr')
parser.add_argument('-long_hmms_x_long_reads_cell','-lh_lr')
parser.add_argument('-short_hmms_x_long_reads_cell','-sh_lr')
parser.add_argument('-long_hmms_x_reads_vir','-lh_rv')
parser.add_argument('-short_hmms_x_reads_vir','-sh_rv')
parser.add_argument('-version', action='store_true')
parser.add_argument('-h', '-help', action='store_true')
arguments = parser.parse_args()

def print_help():
  import os
  mandatory={"-i <hmm_file>":"Profile HMMs to be validated.", "-db_type <'l'|'long'|'s'|'short'>":"Type of database. If -db_type = 'short' or 's', -cell_org_db_frag is mandatory. Else if -db_type = 'long' or 'l', -cell_org_db and -vir_db are mandatory", "-model_type <'l'|'long'|'s'|'short'>":"Type of model.", "-cell_org_db <fasta_file>":"Cellular organisms database, mandatory if -db_type = 'long' or 'l'.","-cell_org_db_frag <fasta_file>":"Cellular organisms fragmented database, mandatory if -db_type = 'short' or 's'.", "-vir_db <directory>":"Directory containing viral sequences, mandatory if -db_type = 'long' or 'l'.", "-conf <file>":"Configuration file (parameters in format 'parameter=value', one parameter per line and parameter names equal to those on the command line)."}
  optional={"-out <name>":"Name of output directory (default: hmm_validated).", "-pt <decimal>":"Maximum percentage ratio (default = 80), must be between 0 and 100. Given similarity results against cell organisms (OrgCellDB) and viral sequences (VirusDB), highest score obtained for OrgCellDB must be lower than pt% of the lowest score obtained against Virus, otherwise the model is not accepted.", "-pd <decimal>":"Minimum accepted percentage detection rate of viral sequences (default = 80) after applying new cutoff scores to the HMMs profile, must be between 0 and 100.", "-sh_sr <string>":"For -model_type = 'short' or 's' and -db_type = 'short' or 's', hmm-prospector cell_org run parameters.","-sh_lr <string>":"For -model_type = 'short' or 's' and -db_type = 'long' or 'l', hmm-prospector cell_org run parameters.", "-lh_sr <string>":"For -model_type = 'long' or 'l' and -db_type = 'short' or 's', hmm-prospector cell_org run parameters.", "-lh_lr <string>":"For -model_type = 'long' or 'l' and -db_type = 'long' or 'l', hmm-prospector cell_org run parameters.", "-sh_rv <string>":"For -model_type = 'short' or 's', hmm-prospector virus run parameters.", "-lh_sr <string>":"For -model_type = 'long' or 'l', hmm-prospector virus run parameters."}
    
  try:
    size = int(str(os.get_terminal_size()).split("columns=")[-1].split(",")[0])
  except:
    size=None
  else:
    pass
  
  if size==None:
    try:
      size = int(subprocess.check_output("tput cols", shell=True))
    except:
      size=None
    else:
      pass
  
  string=''
  
  header = ['(c) 2022. Arthur Gruber & Giuliana Pola']
  header.append('For more information access: https://github.com/GiulianaPola/model_validator')
  header.append("\nUsage:")
  header.append("model_validator.py -conf <configuration_file>")
  header.append("model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'short' -cell_org_db_frag <fasta_file> <optional parameters>")
  header.append("model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'long' -cell_org_db <fasta_file> -vir_db <directory> <optional parameters>")
  header.append("\nMandatory parameters:")
  
  if size==None:
    string=''.join(header)
  else:
    for text in header:
      line=''
      while not text=='':
        i=-1
        if len(text)<=size:
          line+=text
          text=''
        else:
          part=text[0:size]
          i=part.rindex(" ")
          line+=text[0:i]  
          text=text[i+1:]
        string+=line
        string+="\n"
        line=''
  
  keys=list(mandatory.keys())
  keys.extend(list(optional.keys()))
  values=list(mandatory.values())
  values.extend(list(optional.values()))
  maxk=0
  maxv=0
  
  for key in keys:
    if len(key)>maxk:
      maxk=len(key)
  
  for key in sorted(mandatory.keys()):
    line=''
    line+='{message: <{width}} '.format(
    message=key,
    fill=' ',
    align='<',
    width=maxk,
    )
    if size==None:
      line+=mandatory[key]
      line+="\n"
      string+=line
    else:
      rest=size-maxk-1
      text=mandatory[key]
      while not text=='':
        i=-1
        if len(text)<=rest:
          if line=='':
            line+=" "*(maxk+1)
          line+=text
          text=''
        else:
          part=text[0:rest]
          i=part.rindex(" ")
          if line=='':
            line+=" "*(maxk+1)
          line+=text[0:i]  
          text=text[i+1:]
        string+=line
        string+="\n"
        line=''
  
  string += "\nOptional parameters:\n"
  
  for key in sorted(optional.keys()):
    line=''
    line+='{message: <{width}} '.format(
    message=key,
    fill=' ',
    align='<',
    width=maxk,
    )
    if size==None:
      line+=optional[key]
      line+="\n"
      string+=line
    else:
      rest=size-maxk-1
      text=optional[key]
      while not text=='':
        i=-1
        if len(text)<=rest:
          if line=='':
            line+=" "*(maxk+1)
          line+=text
          text=''
        else:
          part=text[0:rest]
          i=part.rindex(" ")
          if line=='':
            line+=" "*(maxk+1)
          line+=text[0:i]  
          text=text[i+1:]
        string+=line
        string+="\n"
        line=''

  print(string)


def rename(i,name,typ):
  path=''
  if '/' in name:
    path=os.path.split(name)[0]
    name=os.path.split(name)[1]
  newname=os.path.join(path, name)
  if typ=='dir':
    while os.path.isdir(newname):
      i+=1
      newname=os.path.join(path, str(name+'_'+str(i)))
  elif typ=='file':
    while os.path.isfile(newname):
      i+=1
      newname=os.path.join(path, str(name+'_'+str(i)))
  return newname

def validate_conf(conf):
    valid=True
    args=dict()
    if not os.path.isfile(conf):
        print("Configuration file (-conf) doesn't exist!")
        valid=False
    else:
        try:
            cfile=open(conf,"r",encoding='utf-8')
        except:
            print("Configuration file (-conf) couldn't be opened!")
            valid=False
        else:
            try:
              text=str(cfile.read()).replace(' = ','=')
            except Exception as e:
              print(e)
              print("Configuration file (-conf) is not encoded in utf-8!")
              valid=False
            else:
              last=""
              for line in sorted(text.split('\n'), reverse=True):
                  last+=line
                  if "=" in last:
                    index=last.find("=")
                    value=last[index+1:]
                    key=last[0:index].lower()
                    if key=='input_file' or key=='i':
                        args['i']=value.replace('"','').strip()
                        last=""
                    elif key=='output' or key=='out' or key=='o':
                        args['out']=value
                        last=""
                    elif key=='db_type' or key=='dbtype':
                        args['db_type']=value
                        last=""
                    elif key=='modeltype' or key=='model_type':
                        args['model_type']=value
                        last=""
                    elif key=='cell_org_db':
                        args['cell_org_db']=value.replace('"','').strip()
                        last=""
                    elif key=='cell_org_db_frag':
                        args['cell_org_db_frag']=value.replace('"','').strip()
                        last=""
                    elif key=='vir_db':
                        args['vir_db']=value.replace('"','').strip()
                        last=""
                    elif key=='vir_db_frag':
                        args['vir_db_frag']=value
                        last=""
                    elif key=='pt':
                        args['pt']=value
                        last=""
                    elif key=='pd':
                        args['pd']=value
                        last=""
                    elif key in ['long_hmms_x_short_reads_cell','lh_sr']:
                        args['long_hmms_x_short_reads_cell']=value
                        last=""
                    elif key in ['short_hmms_x_short_reads_cell','sh_sr']:
                        args['short_hmms_x_short_reads_cell']=value
                        last=""
                    elif key in ['long_hmms_x_long_reads_cell','lh_lr']:
                        args['long_hmms_x_long_reads_cell']=value
                        last=""
                    elif key in ['short_hmms_x_long_reads_cell','sh_lr']:
                        args['short_hmms_x_long_reads_cell']=value
                        last=""
                    elif key in ['long_hmms_x_reads_vir','lh_rv']:
                        args['long_hmms_x_reads_vir']=value
                        last=""
                    elif key in ['short_hmms_x_reads_vir','sh_rv']:
                        args['short_hmms_x_reads_vir']=value
                        last=""
    return valid, args

def validate_args(args):
    valid=True
    if not 'i' in args.keys():
        print("Missing profile HMMs to be validated (-i)!")
        valid=False
    else:
        if os.path.isfile(args['i'].replace('"','')):
            param["i"]=os.path.realpath(args['i'].replace('"',''))
        else:
            print(args['i'].replace('"',''))
            print("Profile HMMs file (-i) doesn't exist!")
            valid=False
    
    if not 'db_type' in args.keys():
        print("Missing database type (-db_type)!")
        valid=False
    elif args['db_type'].lower()=='short' or args['db_type'].lower()=='s':
        args['db_type']="short"
        param['db_type']="short"
        if 'cell_org_db_frag' not in args.keys():
            print("Missing organism cell fragmented database (-cell_org_db_frag)!\n")
            valid=False
        else:
            if not os.path.isfile(args['cell_org_db_frag']):
                print(args['cell_org_db_frag'].replace('"',''))
                print("Organism cell fragmented database (-cell_org_db_frag) doesn't exist!")
                valid=False
            else:
                param["cell_org_db"]=os.path.realpath(args['cell_org_db_frag'])
    elif args['db_type'].lower()=='long' or args['db_type'].lower()=='l':
        param['db_type']="long"
        if 'cell_org_db' not in args.keys():
            print("Missing organism cell database (-cell_org_db)!\n")
            valid=False
        else:
            if not os.path.isfile(args['cell_org_db']):
                print("Organism cell database (-cell_org_db) doesn't exist!")
                valid=False
            else:
                param["cell_org_db"]=os.path.realpath(args['cell_org_db'])
        if 'vir_db' not in args.keys():
            print("Missing virus database (-vir_db)!\n")
            valid=False
        else:
            if not os.path.isdir(args['vir_db']):
                print("Virus database (-vir_db) doesn't exist!")
                valid=False
            else:
                param["vir_db"]=os.path.realpath(args['vir_db'])
    else:
        print("Database type (-db_type) must be 'short' or 'long'!")
        valid=False

    if not 'model_type' in args.keys():
        print("Missing profile HMMs type (-model_type)!")
        valid=False
    elif args['model_type'].lower() in ['s','short']:
        param['model_type']="short"
    elif args['model_type'].lower() in ['l','long']:
        param['model_type']="long"
    else:
        print("Profile HMMs type (-model_type) must be 'short' or 'long'!")
        valid=False

    if valid==True:
        if 'reads_cell' in ','.join(list(args.keys())):
          if "', " in str(args).split("reads_cell': ")[1]:
            param['reads_cell']=str(args).split("reads_cell': ")[1].split(",")[0]
          else:
            param['reads_cell']=str(args).split("reads_cell': ")[1].split("}")[0]
        if 'reads_vir' in ','.join(list(args.keys())):
          if "', " in str(args).split("reads_cell': ")[1]:
            param['reads_vir']=str(args).split("reads_vir': ")[1].split(",")[0]
          else:
            param['reads_vir']=str(args).split("reads_vir': ")[1].split("}")[0]
        if param['db_type']=='long':
          if 'reads_cell' in param.keys() and 'reads_vir' not in param.keys():
            param['reads_vir']=param['reads_cell']
          elif 'reads_vir' in param.keys() and 'reads_cell' not in param.keys():
            param['reads_cell']=param['reads_vir']
          elif 'reads_cell' not in param.keys() and 'reads_vir' not in param.keys():
            if param['model_type']=='short':
              args['short_hmms_x_long_reads_cell']="-cpu 20 -n 10 -r no -sv 1.0"
              args['short_hmms_x_reads_vir']="-cpu 20 -n 10 -r no -sv 1.0"
              param['reads_cell']="-cpu 20 -n 10 -r no -sv 1.0" 
              param['reads_vir']="-cpu 20 -n 10 -r no -sv 1.0"
            elif param['model_type']=='long':
              args['long_hmms_x_long_reads_cell']="-cpu 20 -n 10 -r no -sv 0.6" 
              args['long_hmms_x_reads_vir']="-cpu 20 -n 10 -r no -sv 0.6" 
              param['reads_cell']="-cpu 20 -n 10 -r no -sv 0.6"
              param['reads_vir']="-cpu 20 -n 10 -r no -sv 0.6"
        elif param['db_type']=='short':
          if 'reads_cell' not in param.keys():
            if param['model_type']=='short':
              args['short_hmms_x_short_reads_cell']="-cpu 20 -n 10 -rl c -sv 1.5"
              param['reads_cell']="-cpu 20 -n 10 -rl c -sv 1.5"
            elif param['model_type']=='long':
              args['long_hmms_x_short_reads_cell']="-cpu 20 -n 10 -rl c -sv 0.8"
              param['reads_cell']="-cpu 20 -n 10 -rl c -sv 0.8"

        if 'pt' in args.keys():
          if ',' in args['pt']:
              pt=args['pt'].replace(',','.')
          else:
              pt=args['pt']
          try:
              pt=float(pt)
          except Exception as e:
              print(e)
              print("Maximum percentage ratio (-pt) is not a decimal!")
              valid=False
          else:
              if pt>0 and pt<=100:
                  param['pt']=pt
              else:
                  print("Maximum percentage ratio (-pt) must be between 0 and 100!")
                  valid=False
        else:
          param['pt']=80
          args['pt']=80

        if 'pd' in args.keys():
          if ',' in args['pd']:
              pd=args['pd'].replace(',','.')
          else:
              pd=args['pd']
          try:
              pd=float(pd)
          except Exception as e:
              print(e)
              print("Minimum detection rate (-pd) is not a decimal!")
              valid=False
          else:
              if pd>0 and pd<=100:
                  param['pd']=pd
              else:
                  print("Minimum detection rate (-pd) must be between 0 and 100!")
                  valid=False
        else:
          args['pd']=80
          param['pd']=80

        if not 'out' in args.keys():
          out=os.path.join(call,'hmm_validated')
        else:
          head_tail = os.path.split(args['out'])
          if os.path.exists(head_tail[0]):
            out=args['out']
          else:
            out=os.path.join(call,head_tail[1])
        if os.path.isdir(out):
          print("Output directory '{}' already exist!".format(os.path.split(out)[-1]))
          out=rename(0,out,'dir')
        try:
          os.mkdir(out)
        except:
          print("Output directory '{}' couldn't be created!".format(os.path.split(out)[-1]))
          valid=False
        else:
          print("Creating output directory '{}'...".format(os.path.split(out)[-1]))
          param['out']=out
          args['out']=out
    return valid,param

def cmd_prospector(parameters,db,i,o,name):
  #log.write("CMD_PROSPECTOR\nparameters={}\ndb={}\ni={}\no={}\nname={}\n".format(parameters,db,i,o,name))
  if ", " in parameters:
    parameters=parameters.split(", ")[0]
  parameters=" "+parameters.strip("'").strip('"').strip("'")
  cmd="nohup hmm-prospector.pl"
  cmd+=" -d {}".format(db)
  cmd+=" -i {}".format(i)
  if "hmm-prospector.pl" in parameters:
    parameters=parameters.split('hmm-prospector.pl')[1]
  for arg in parameters.split(' -'):
    if not arg.strip()=="":
     key=arg.split()[0]
     value=arg.split()[1]
     if not key=="o":
       if " -{} ".format(key) not in cmd:
         cmd+=" -{} {}".format(key, value)
  if name=='cell_org':
    cmd+=" -o hmm-prospector"
  else:
    cmd+=" -o {}".format(o)
  cmd+=" 1>{}_prosp_out 2>{}_prosp_er &".format(os.path.split(o)[-1],os.path.split(o)[-1])
  if name=="cell_org":
    log.write("hmm-prospector for cell_org: {}\n".format(cmd))
  else:
    log.write("'{}' hmm-prospector for {}: {}\n".format(os.path.split(o)[-1],name,cmd))
  return cmd

def run_prospector(db,i,o,parameters,name,look,cmd):
  positive=[]
  valid=True
  os.system(cmd)
  while not os.path.isfile(look):
    pass
  output=open(look,"r")
  process=os.popen("ps -u {} --sort +start_time".format(os.getlogin())).read()
  table="{}/hmm-prospector/table2.csv".format(o)
  while "hmm-prospector." in process and "Done." not in output.read():
    output.close()
    output=open(look,"r")
    process=os.popen("ps -u {} --sort +start_time".format(os.getlogin())).read()
  output.close()
  return valid

def parse_table(name,o,i,match):
  valid=True
  positive=[]
  table=os.path.join(o,name,"hmm-prospector","table2.csv")
  log.write("Match table: {}\n".format(table)) 
  if not os.path.isfile(table):
    print("Match {} table doesn't exist!".format(name))
    log.write("Match {} table doesn't exist!\n".format(name))
    valid=False
  else:
   try:
     tfile=open(table,'r')
   except Exception as e: 
     log.write("{}\n".write(e))
     print("Match {} table couldn't be opened!".format(name))
     log.write("Match {} table couldn't be opened!\n".format(name))
     valid=False
   else:
     for line in tfile.readlines():
       hmm=line.split("\t")[0]
       if not hmm=="HMM":
         positive.append(hmm)
     # cut -f1 -d$'\t' "/home/geninfo/gpola/model_validator/short_dbXshort_model/cell_org/hmm-prospector/table2.csv" > match
     tfile.close()
     log.write("{} matched HMMs for {}.\n".format(len(positive),name))
  return valid,positive
  
def get_names(file):
  names=[]
  lines=str(subprocess.check_output("grep ' {} ' {}".format(name,file), shell=True)).split("\\n")
  if "'" in lines:
    lines.remove("'")
  for line in lines:
    columns=line.split(" ")
    while '' in columns:
      columns.remove('')
    while '-' in columns:
      columns.remove('-')
    if len(columns)>=2:
      names.append(columns[1])
  return names
     
def get_diff(hmm,positive,match,name):
  if match=="Valid":
    mismatch="Invalid"
  elif match=="Invalid":
    mismatch="Valid"
  elif match=="Match":
    mismatch="No match"
  elif match=="No match":
    mismatch="Match"
  valid=True
  negative=[]
  try:
     ifile=open(hmm,"r")
  except Exception as e: 
     log.write("{}\n".write(e))
     valid=False
     log.write("Profile HMMs to be validated (-i) couldn't be opened!\n")
     print("Profile HMMs to be validated (-i) couldn't be opened!")
  else:
     text="\n".join(sorted(ifile.readlines()))
     text='\n'.join(list(text.split("\nNAME  ")[1:-1]))+'\n'+str(text.split("\nNAME  ")[-1].split("\n")[0])
     negative=text.replace("\n\n","\n").split("\n")
     log.write("Total number of HMMs: {}\n".format(len(negative)))
     if not positive==[]:
       for model in positive:
         if model in negative:
           negative.remove(model)
     if "" in negative:
       negative.remove("")
#  if not positive==[]:
#    if not name=="":
#      log.write("{} {} HMMs for {}: {}\n".format(len(positive),match.lower(),name,positive))
#    else:
#      log.write("{} {} HMMs: {}\n".format(len(positive),match.lower(),positive))
#  if not negative==[]:   
#    if not name=="":
#      log.write("{} {} HMMs for {}: {}\n".format(len(negative),mismatch.lower(),name,negative))
#    else:
#      log.write("{} {} HMMs: {}\n".format(len(negative),mismatch.lower(),negative))
  return valid,negative

def list_HMMs(name,match,o,positive):
  valid=True
  listHMMs=''
  if positive==[]:
    if not name=='':
      print("List of {} HMMs for {} is empty!".format(match.lower(),name))
      log.write("List of {} HMMs for {} is empty!\n".format(match.lower(),name))
    else:
      print("List of {} HMMs is empty!".format(match.lower()))
      log.write("List of {} HMMs is empty!\n".format(match.lower()))
  else:
    try:
      if match=="Match":
        listHMMs=os.path.join(o,name,"match_HMMs_list.txt")
      elif match=="Valid" or match=="Invalid":
        if not os.path.isdir(os.path.join(o,"results")):
          os.mkdir(os.path.join(o,"results"))
        listHMMs=os.path.join(o,"results","{}_HMMs_list.txt".format(match.lower()))
      mfile=open(listHMMs,"w")
      mfile.write('\n'.join(positive))
      mfile.close()
    except Exception as e: 
      log.write("{}\n".format(e))
      if not name=='':
        print("List of {} HMMs for {} could not be written!".format(match.lower(),name))
        log.write("List of {} HMMs for {} could not be written!\n".format(match.lower(),name))
      else:
        print("List of {} HMMs could not be written!".format(match.lower()))
        log.write("List of {} HMMs could not be written!\n".format(match.lower()))
      valid=False
    else:
      pass
  return valid,listHMMs
  
def validate_database(d):
  valid=True
  database=''
  dtext=''
  try:
    dfile=open(d,"r")
  except Exception as e: 
    log.write("{}\n".format(e))
    print("Profile HMM dataset couldn't be opened!")
    valid=False
  else:
    dtext=dfile.read()
    if dtext.isspace() or len(dtext)==0:
      print("Profile HMM dataset is empty!")
      valid=False
    else:
      database=dtext
  return valid,database

def fetch_hmm(namelist,d,o):
  valid=True
  if not os.path.isdir(o):
    try:
      os.mkdir(o)
    except Exception as e: 
      log.write("{}\n".format(e))
      print("'{}' directory was not created!".format(os.path.split(o)[-1]))
      log.write("'{}' directory was not created!\n".format(os.path.split(o)[-1]))
      valid=False
    else:
      print("'{}' directory was created!".format(os.path.split(o)[-1]))
      log.write("'{}' directory was created!\n".format(os.path.split(o)[-1])) 
  hmms=[]
  if valid==True:
    #valid,namelist=validate_list(i)
    if valid==True:
      valid,database=validate_database(d)
      if valid==True:
        missing=[]
        selected=[]
        for name in namelist:
          if name not in database:
            missing.append(name)
            log.write("{} was not found in HMM file (-i)!".format(name))
          else:
            selected.append(name)
            head=database.split("NAME  "+name)[0]
            tail=database.split("NAME  "+name)[1]
            hmm="HMMER"+head.split("HMMER")[-1]+"NAME  "+name+tail.split("//")[0]+"//"
            hmms.append(hmm)
            try:
              hmmfile=open(os.path.join(o, '{}.hmm'.format(name)),'w')
              hmmfile.write("")
            except Exception as e: 
              log.write("{}\n".format(e))
              print("'{}' HMM file was not created!".format(name))
              log.write("'{}' HMM file was not created!\n".format(name))
              valid=False
            else:
              hmmfile.write(hmm)
              hmmfile.close()
        if not missing==[]:
          log.write("Missing HMMs: {}\n".format(missing))
        if selected==[]:
          print("All selected HMMs are missing!")
          log.write("All selected HMMs are missing!\n")
          valid=False
        else:
          try:
            filename=os.path.basename(os.path.normpath(o))
            selectedfile=open(os.path.join(o, '{}.hmm'.format(filename)),'w')
            selectedfile.write("")
          except Exception as e: 
            log.write("{}\n".format(e))
            print("Selected HMM file was not created!")
            log.write("Selected HMM file was not created!\n")
            valid=False
          else:
            selectedfile.write('\n'.join(hmms))
  return valid

def call_fetch(name,o,i,match,positive,both):
  if match[0]=="Valid":
    mismatch="Invalid"
  elif match[0]=="Invalid":
    mismatch="Valid"
  elif match[0]=="Match":
    mismatch="No match"
  elif match[0]=="No match":
    mismatch="Match"
  valid=True
  negative=[]  
  valid,listHMMs=list_HMMs(name,match[0],o,positive)
  if valid and not listHMMs=='':
    if match[0]=="Match":
      out=os.path.join(o,name,"models")
    elif match[0]=="Valid" or match[0]=="Invalid":
      if not os.path.isdir(os.path.join(o,"results")):
        os.mkdir(os.path.join(o,"results"))
      out=os.path.join(o,"results",match[0].lower())
    start=datetime.now()
    positive=sorted(list(set(positive)))
    if not fetch_hmm(positive,i,out):
      error.extend(positive)
      valid=False
      if match[0] in ["Valid","Invalid"]:
        print("ERROR: {} HMMs couldn't be fetched!".format(match[0].capitalize()))
        log.write("ERROR: {} HMMs couldn't be fetched!\n".format(match[0].capitalize()))
      else:
        print("ERROR: HMMs matching {} couldn't be fetched!".format(name))
        log.write("ERROR: HMMs matching {} couldn't be fetched!\n".format(name))
    else:
      if match[0] =="Match":
        log.write("Runtime to fetch the HMMs matching {}: {}\n".format(name,datetime.now()-start))
        print("{} HMMs matching {} were fetched!".format(len(positive),name))
        log.write("{} HMMs matching {} were fetched!\n".format(len(positive),name))
        log.write("{} {} match HMMs: {}\n".format(len(positive),name,positive))
      elif match[0]=="Invalid" or match[0]=="Valid":
        log.write("{} {} HMMs: {}\n".format(len(positive),match[0].capitalize(),positive))
        print("Runtime to fetch {} HMMs: {}".format(match[0].lower(),datetime.now()-start))
        log.write("Runtime to fetch {} HMMs: {}\n".format(match[0].lower(),datetime.now()-start))
        print("{} {} HMMs were fetched!".format(len(positive),match[0].lower()))
        log.write("{} {} HMMs were fetched!\n".format(len(positive),match[0].lower()))
  if both==True:
    valid,negative=get_diff(i,positive,match[0],name)
    if negative==[]:
      if not match[1]==mismatch:
        if not name=="":
          log.write("{} HMMs list for {} is empty!\n".format(mismatch.capitalize(),name))
          print("{} HMMs list for {} is empty!".format(mismatch.capitalize(),name))
        else:
          log.write("{} HMMs list is empty!\n".format(mismatch.capitalize()))
          print("{} HMMs list is empty!".format(mismatch.capitalize()))
      else:
        log.write("{} HMMs list is empty!\n".format(match[1].capitalize()))
        print("{} HMMs list is empty!".format(match[1].capitalize()))
    else:
        negative=sorted(list(set(negative)))
        if not os.path.isdir(os.path.join(o,"results")):
          try:
            os.mkdir(os.path.join(o,"results"))
          except:
            error.extend(negative)
            print("ERROR: 'results' directory was not created!")
            log.write("ERROR: 'results' directory was not created!\n")
            valid=False
          else:
            print("'results' directory was created!")
            log.write("'results' directory was created!\n")
        #log.write("{}:{}\n".format(mismatch,negative))
        valid,listHMMs=list_HMMs("",match[1],o,negative)
        if valid and not listHMMs=='':
          if match[1]=="Valid" or match[1]=="Invalid":
            out=os.path.join(o,"results",match[1].lower())
            start=datetime.now()
            if not fetch_hmm(negative,i,out):
              error.extend(negative)
              valid=False
              if not match[1]==mismatch:
                if not name=="":
                  print("ERROR: {} HMMs for {} couldn't be fetched!".format(mismatch.capitalize(),name))
                  log.write("ERROR: {} HMMs for {} couldn't be fetched!\n".format(mismatch.capitalize(),name))
                else:
                  print("ERROR: {} HMMs couldn't be fetched!".format(mismatch.capitalize()))
                  log.write("ERROR: {} HMMs couldn't be fetched!\n".format(mismatch.capitalize()))
              else:
                print("ERROR: {} HMMs couldn't be fetched!".format(match[1].capitalize()))
                log.write("ERROR: {} HMMs couldn't be fetched!\n".format(match[1].capitalize()))
            else:
              if not match[1]==mismatch:
                if not name=="":
                  print("{} {} HMMs for {} were fetched as {}!".format(len(negative),mismatch.lower(),name,match[1].lower()))
                  log.write("{} {} HMMs for {} were fetched as {}!\n".format(len(negative),mismatch.lower(),name,match[1].lower()))
                  log.write("{} for {}:{}\n".format(mismatch.capitalize(),name,negative))
                else:
                  print("{} {} HMMs were fetched as {}!".format(len(negative),mismatch.lower(),match[1].lower()))
                  log.write("{} {} HMMs were fetched as {}!\n".format(len(negative),mismatch.lower(),match[1].lower()))
                  log.write("{}:{}\n".format(mismatch.capitalize(),negative))
              else:
                print("{} {} HMMs were fetched!".format(len(negative),match[1].lower()))
                log.write("{} {} HMMs were fetched!\n".format(len(negative),match[1].lower()))
                log.write("{}:{}\n".format(match[1].capitalize(),negative))
  return valid
  
def get_path(path,name):
  import os
  folders=os.listdir(path)
  try:
    idx = [item.lower() for item in folders].index(name.lower())
  except Exception as e: 
    #log.write("{}\n".format(e))
    return False
  else:
    return os.path.join(path,folders[idx])
  
def get_vir_db(table,folder,models):
  valid=[]
  error=[]
  taxons=[]
  vir_db=dict()
  fam_dir=''
  prot_dir=''
  headers=str(subprocess.check_output("head {} -n 1".format(table), shell=True))[2:-1].replace("\\n","").split("\\t")
  taxons=str(subprocess.check_output("cut -f{} {} | sort -u".format(headers.index("Taxon")+1,table), shell=True))[2:-1].split("\\n")
  names=str(subprocess.check_output("grep 'vHMM_' {} | cut -f{} | sort -u".format(table,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
  taxons.remove("")
  taxons.remove("Taxon")
  for taxon in taxons:
    family=None
    names=str(subprocess.check_output("grep '{}' {} | cut -f{} | sort -u".format(taxon,table,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
    names.remove("")
    family=''
    rank=''
    txid=str(subprocess.check_output('grep "TXID" "{}/{}.hmm"'.format(models,names[0]), shell=True)).split("\\t")[-1].split("\\n")[0]
    Entrez.email = "giulianapola@usp.br"
    handle=None
    for n in range(3):
      try:
        handle = Entrez.efetch(db="Taxonomy", id=txid, retmode="xml")
        break
      except:
        handle=None
        time.sleep(1)
    if handle==None:
      print("ERROR: Failure in the connection with NCBI!".format(taxon,",".join(names)))
      for name in names:
        log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
        log.write("'{}' taxon: {}\n".format(name,taxon))
        log.write("'{}' ERROR: Failure in the connection with NCBI!\n".format(name))
      error.extend(names)
    else:
      records = Entrez.read(handle)[0]
      if not "Rank" in records.keys():
        error.extend(names)
        print("ERROR: Taxon '{}' rank could not be found: {}".format(taxon,",".join(names)))
        for name in names:
          log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
          log.write("'{}' taxon: {}\n".format(name,taxon))
          log.write("'{}' ERROR: Taxon '{}' rank could not be found!\n".format(name,taxon))
      else:
        rank=records["Rank"]
        handle.close()
        if rank=='family':
          family=taxon
        elif not rank in ['genus','subfamily']:
          error.extend(names)
          print("ERROR: Taxon '{}' rank is not genus or family or subfamily: {}".format(taxon,",".join(names)))
          for name in names:
            pass
            log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
            log.write("'{}' taxon: {}\n".format(name,taxon))
            log.write("'{}' taxon rank: {}\n".format(name,rank))
            log.write("'{}' ERROR: Taxon '{}' rank is not genus or family or subfamily!\n".format(name,taxon))
        else:
          if "LineageEx" in records.keys():
            for item in reversed(records["LineageEx"]):
               if item['Rank']=='family':
                 family=item['ScientificName']
        if family==None:
          error.extend(names)
          print("ERROR: Taxon '{}' family was not found: {}".format(taxon,",".join(names)))
          for name in names:
            log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
            log.write("'{}' taxon: {}\n".format(name,taxon))
            log.write("'{}' taxon rank: {}\n".format(name,rank))
            log.write("'{}' ERROR: Taxon '{}' family was not found!\n".format(name,taxon))
        else:
          fam_dir=get_path(folder,family)
          if fam_dir==False:
              names=str(subprocess.check_output("grep '{}' {} | cut -f{} | sort -u".format(taxon,table,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
              names.remove("")
              error.extend(names)
              print("ERROR: Missing family '{}' folder: {}".format(family,",".join(names)))
              for name in names:
                pass
                log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
                log.write("'{}' taxon: {}\n".format(name,taxon))
                log.write("'{}' taxon rank: {}\n".format(name,rank))
                log.write("'{}' family: {}\n".format(name,family))
                log.write("'{}' ERROR: Missing family folder: {}\n".format(name,os.path.join(folder,family)))
          else:
            proteins=str(subprocess.check_output("grep '{}' {} | cut -f4 | sort -u".format(taxon,table), shell=True))[2:-1].split("\\n")
            proteins.remove("")
            for protein in proteins:
              prot_file=''
              prot_dir=get_path(fam_dir,protein)
              if prot_dir==False:
                n=0
                for filename in os.listdir(fam_dir):
                  if os.path.isdir(os.path.join(fam_dir,filename)) and filename.replace(" ","_").lower() in protein.replace(" ","_").lower():
                    if len(filename)>n:
                      prot_dir=os.path.join(fam_dir,filename)
                      n=len(filename)
                  if os.path.isdir(os.path.join(fam_dir,filename)) and protein.replace(" ","_").lower() in filename.replace(" ","_").lower():
                    if len(protein)>n:
                      prot_dir=os.path.join(fam_dir,filename)
                      n=len(protein)
              if prot_dir==False:
                names=str(subprocess.check_output("grep '{}' {} | grep '{}' | cut -f{} | sort -u".format(taxon,table,protein,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
                names.remove("")
                error.extend(names)
                print("ERROR: Missing protein folder '{}/{}': {}".format(os.path.split(fam_dir)[-1],protein.replace(" ","_").lower(),",".join(names)))
                for name in names:
                  pass
                  log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
                  log.write("'{}' taxon: {}\n".format(name,taxon))
                  log.write("'{}' taxon rank: {}\n".format(name,rank))
                  log.write("'{}' family: {}\n".format(name,family))
                  log.write("'{}' protein: {}\n".format(name,protein))
                  log.write("'{}' ERROR: Missing protein folder: {}\n".format(name,os.path.join(fam_dir,protein.replace(" ","_").lower())))
              else:
                seq_dir=None
                if rank=='subfamily':
                  path=os.path.join("full-length","tabajara_narrow_subfamilies","fastas")
                  if os.path.isdir(os.path.join(prot_dir,"full-length","tabajara_narrow_subfamilies","fastas")):
                    seq_dir=os.path.join(prot_dir,"full-length","tabajara_narrow_subfamilies","fastas")
                else:
                  path="fasta_sequences"
                  if os.path.isdir(os.path.join(prot_dir,"fasta_sequences")):
                    seq_dir=os.path.join(prot_dir,"fasta_sequences")
                if seq_dir==None:
                  names=str(subprocess.check_output("grep '{}' {} | grep '{}' | cut -f{} | sort -u".format(taxon,table,protein,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
                  names.remove("")
                  error.extend(names)
                  print("ERROR: Missing sequences folder '{}': {}".format(path,",".join(names)))
                  for name in names:
                    pass
                    log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
                    log.write("'{}' taxon: {}\n".format(name,taxon))
                    log.write("'{}' taxon rank: {}\n".format(name,rank))
                    log.write("'{}' family: {}\n".format(name,family))
                    log.write("'{}' protein: {}\n".format(name,protein))
                    log.write("'{}' ERROR: Missing sequences folder: {}\n".format(name,os.path.join(prot_dir,path)))
                else:
                  names=str(subprocess.check_output("grep '{}' {} | grep '{}' | cut -f{} | sort -u".format(taxon,table,protein,headers.index("HMM")+1), shell=True))[2:-1].split("\\n")
                  names.remove("")
                  if rank=='family':
                    path=os.path.join(seq_dir,"{}_no_redundancy.fasta".format(taxon))
                    prot_file=get_path(seq_dir,"{}_no_redundancy.fasta".format(taxon))
                  else:
                    path=os.path.join(seq_dir,"{}.fasta".format(taxon))
                    prot_file=get_path(seq_dir,"{}.fasta".format(taxon))
                  if prot_file==False:
                    error.extend(names)
                    path_list=os.path.split(path)
                    i=path_list.index(family)
                    print("ERROR: Missing protein file '{}': {}".format(os.path.split(prot_file)[-1],",".join(names)))
                    for name in names:
                      pass
                      log.write("\n{}\n'{}' match cell_org!\n".format(name,name))
                      log.write("'{}' taxon: {}\n".format(name,taxon))
                      log.write("'{}' taxon rank: {}\n".format(name,rank))
                      log.write("'{}' family: {}\n".format(name,family))
                      log.write("'{}' protein: {}\n".format(name,protein))
                      log.write("'{}' ERROR: Missing protein file: {}\n".format(name,os.path.join(path_list[i:])))
                  else:
                    for name in names:
                      valid.append(name)
                      vir_db[name]=prot_file
                      taxonomy[name]={'Rank':rank,'Taxon':taxon,'Protein':protein}
                      if not family=='':
                        taxonomy[name]['Family']=family
  error=sorted(list(set(error)))
  if len(error)>0:
    print("ERROR: {} HMM models had missing taxon folder or missing protein file!".format(len(error)))
    log.write("\nERROR: Some HMM models were invalidated because of missing taxon folder or missing protein file!\n")
    log.write("{} HMMs with missing files: {}\n".format(len(error),",".join(error)))
    log.write("{} HMMs with files found: {}\n".format(len(valid),",".join(valid)))
  else:
    print("All taxon folders and protein files from {} HMM models were found!".format(len(valid)))
    log.write("\nAll taxon folders and protein files were found!\n")
    log.write("{} with files found: {}\n".format(len(valid),",".join(valid)))
  if len(vir_db)==0:
    return False,error
  else:
    return vir_db,error

def get_score(name,file,mode,index):
  if mode=="min":
    complement="head -n1"
  elif mode=="max":
    complement="tail -n1"
  score=float(subprocess.check_output("grep '{}' {} | cut -f{} | sort -n | {}".format(name,file,index,complement), shell=True))
  return score

def compare_scores(cell_org,vir,pt,i,o,name):
  valid=True
  import subprocess
  headers = str(subprocess.check_output("head {} -n 1".format(vir), shell=True)).replace("b'","").replace("\\n","").split("\\t")
  indexes=[]
  indexes.append(headers.index("query_pHMM")+1)
  indexes.append(headers.index("Score")+1)
  vir_hmm=open(i,"r").read()
  new_hmm=''
  oldcutoff=float(vir_hmm.split("CUTOFF SCORE")[-1].split("\n")[0].strip())
  try:
    V = get_score(name,vir,"min",indexes[1])
  except Exception as e: 
    log.write("INVALID: No '{}' vir matches that satisfy reporting thresholds!\n".format(name))
    ##print("INVALID: No '{}' vir matches that satisfy reporting thresholds!".format(name))
    key="invalidated because had no vir matches"
    if not key in finalresults.keys():
      finalresults[key]=[name]
    elif not name in finalresults[key]:
      finalresults[key].append(name)
    valid=False
  else:
    if V<oldcutoff:
      print("ERROR: '{}' min vir score < old cutoff score!".format(name))
      log.write("ERROR: '{}' min vir score < old cutoff score!\n".format(name))
      error.append(name)
    else:
      try:
        C = get_score(name,cell_org,"max",indexes[1])
      except Exception as e: 
        log.write("INVALID: No '{}' cell_org matches that satisfy reporting thresholds!\n".format(name))
        ##print("INVALID: No '{}' cell_org matches that satisfy reporting thresholds!".format(name))
        key="invalidated because had no cell_org matches"
        if not key in finalresults.keys():
          finalresults[key]=[name]
        elif not name in finalresults[key]:
          finalresults[key].append(name)
        invalidated.append(name)
        valid=False
      else:
        if C<oldcutoff:
          print("ERROR: '{}' max cell_org score < old cutoff score!".format(name))
          log.write("ERROR: '{}' max cell_org score < old cutoff score!\n".format(name))
          error.append(name)
        else:
          P = float(V) * pt/100
          if C>P:
            csv.write("{}\t{}\t{}\t{}\tFalse\t\t\t\t\n".format(name,V,C,P))
            log.write("'{}' {:.1f}% (-pt) of min vir score < max cell_org score!\n".format(name,pt))
            log.write("query_pHMM={}\tmin_vir_score={:.1f}\tmax_cell_org_score={:.1f}\tpt%_min_vir_score={:.1f}\n".format(name,float(V),float(C),float(P)))
            newcutoff = 1.2*C
          else:
            newcutoff = (pt/100 * (V - C)) + C
            log.write("'{}' {:.1f}% (-pt) of min vir score > max cell_org score!\n".format(name,pt))
            log.write("query_pHMM={}\tmin_vir_score={:.1f}\tmax_cell_org_score={:.1f}\tpt%_min_vir_score={:.1f}\n".format(name,float(V),float(C),float(P)))
          if vir_hmm=='':
            log.write("ERROR: Profile HMM file '{}.hmm' in cell_org is empty!\n".format(name))
            print("ERROR: Profile HMM file '{}.hmm' in cell_org is empty!".format(name))
            error.append(name)
            valid=False
          else:
            log.write("query_pHMM={}\told_cutoff_score={:.1f}\tnew_cutoff_score={:.1f}\n".format(name,oldcutoff,newcutoff))
            if newcutoff<oldcutoff:
              csv.write("{}\t{}\t{}\t{}\tTrue\t{}\t{}\t{}\tFalse\n".format(name,V,C,P,V-C,newcutoff,oldcutoff))
              print("ERROR: '{}' new cutoff score < old cutoff score!".format(name))
              log.write("ERROR: '{}' new cutoff score < old cutoff score!\n".format(name))
              error.append(name)
              valid=False
            else:
              csv.write("{}\t{}\t{}\t{}\tTrue\t{}\t{}\t{}\tTrue\n".format(name,V,C,P,V-C,newcutoff,oldcutoff))
              new_hmm+=vir_hmm.split("CUTOFF SCORE")[0]+"CUTOFF SCORE    "+"{:.1f}".format(newcutoff)+"\nTAXON"+vir_hmm.split("TAXON")[-1]
              try:
                out=open(o,"w")
                out.write("")
              except:
                log.write("ERROR: Profile HMM file '{}.hmm' with new cutoff score was not created!\n".format(name))
                print("ERROR: Profile HMM file '{}.hmm' with new cutoff score was not created!".format(name))
                error.append(name)
                valid=False
              else:
                out.write(new_hmm)
                out.close()
  return valid
  
def validate_recall(path,name,fasta,pd):
  valid=True
  match=int(str(subprocess.check_output("wc -l {}/recall/{}/table1.csv".format(path,name), shell=True)).replace("b'","").split()[0])-1
  if match<=0:
    #print("INVALID: '{}' recall didn't have any matches!")
    log.write("INVALID: '{}' recall didn't have any matches!\n")
    key="invalidated because recall < (pd%)"
    if not key in finalresults.keys():
      finalresults[key]=[name]
    elif not name in finalresults[key]:
      finalresults[key].append(name)
    valid=False
  else:
    total=int(str(subprocess.check_output("grep '>' {}| wc -l".format(fasta), shell=True)).replace("b'","").replace("\\n'","").split()[0])
    if match/total*100<pd:
      #print("INVALID: '{}' recall detection percentage {:.1f}% <= {:.1f}% (-pd)!".format(name,match/total*100,pd))
      log.write("INVALID: '{}' recall detection percentage {:.1f}% <= {:.1f}% (-pd)!\n".format(name,match/total*100,pd))
      key="invalidated because recall < (pd%)"
      if not key in finalresults.keys():
        finalresults[key]=[name]
      elif not name in finalresults[key]:
        finalresults[key].append(name)
      valid=False
    else:
      log.write("'{}' recall detection percentage {:.1f}% > {:.1f}% (-pd)!\n".format(name,match/total*100,pd))
      key="validated with adjusted cutoff score"
      if not key in finalresults.keys():
        finalresults[key]=[name]
      elif not name in finalresults[key]:
        finalresults[key].append(name)
  return valid 

if not len(sys.argv)>1:
    print_help()
elif arguments.h == True:
    print_help()
elif arguments.version == True:
    print(version) 
else:
    valid=True
    if arguments.conf==None:
        args={k: v for k, v in vars(arguments).items() if v is not None and v is not True and v is not False}     
    else:
        valid,args=validate_conf(arguments.conf)
    if valid:
      valid,param=validate_args(args)
      if valid==False:
        print("ERROR: Input parameters are invalid!")
      if valid==True:
        print("Input parameters are valid!")
        os.chdir(param['out'])
        try:
          log=open(os.path.join(param["out"], 'file.log'),'w')
          log.write('model_validator v{}\n'.format(version))
        except Exception as e: 
          print(e)
          print('ERROR: Log file was not created!')
        else:
          log.write('\nDatetime:\n{}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
          
          log.write('\nWorking directory:\n{}\n'.format(call))
          
          log.write('\nCommand line:\n{}\n'.format(' '.join(sys.argv)))
          
          user=""
          try:
            user=os.getlogin()
          except:
            try:
              user=os.environ['LOGNAME']
            except:
              try:
                user=os.environ['USER']
              except:
                pass
              else:
                pass
            else:
              pass
          else:
            pass
          if not user=="":
            log.write('\nUser:\n{}\n'.format(user))
          
          log.write('\nParameters:\n')
          for key in sorted(args.keys()):
            if not args[key] in [None, False, True]:
              if key in param.keys():
                log.write("{}={}\n".format(key,param[key]))
              else:
                log.write("{}={}\n".format(key,args[key]))
          log.write('\nDataset analysis:\n')
          log.close()
          log=open(os.path.join(param["out"], 'file.log'),'a')
          try:
            if not os.path.isdir(os.path.join(param['out'],'cell_org')):
              os.mkdir(os.path.join(param['out'],'cell_org'))
            os.chdir(os.path.join(param['out'],'cell_org'))
          except Exception as e: 
            print(e)
            print("ERROR: cell_org directory wasn't created!")
            log.write("\nERROR: cell_org directory wasn't created!\n")
          else:
            total=int(str(subprocess.check_output("grep 'NAME' {}| wc -l".format(param['i']), shell=True)).replace("b'","").replace("\\n'","").split()[0])
            print("{} tested profile HMMs".format(total))
            db=param['cell_org_db']
            i=param['i']
            o=os.path.join(param['out'],"cell_org")
            name='cell_org'
            look=os.path.join(param['out'],'cell_org','cell_org_prosp_out')
            start=datetime.now()
            cmd=cmd_prospector(param['reads_cell'],db,i,o,name)
            log.close()
            log=open(os.path.join(param["out"], 'file.log'),'a')
            if run_prospector(db,i,o,param['reads_cell'],name,look,cmd):
              #print("hmm-prospector runtime for {}: {}".format("cell_org",datetime.now()-start))
              log.write("hmm-prospector runtime for {}: {}\n".format("cell_org",datetime.now()-start))
              log.close()
              log=open(os.path.join(param["out"], 'file.log'),'a')
              if param['db_type']=='short':
                try:
                  if not os.path.isdir(os.path.join(param['out'],'results')):
                    os.mkdir(os.path.join(param['out'],'results'))
                  os.chdir(os.path.join(param['out'],'results'))
                except Exception as e: 
                  log.write(e)
                  print("\nERROR: Results directory wasn't created!")
                  log.write("\nERROR: Results directory wasn't created!\n")
                else:
                  print("Results directory was created!")
                  log.write("\nResults directory was created!\n")
                  name='cell_org'
                  o=param['out']
                  i=param['i']
                  match="Invalid"
                  positive=[]
                  valid,positive=parse_table(name,o,i,match)
                  finalresults["positive for OrgCell"]=positive
                  print("{} positive for OrgCell = invalidated".format(len(set(finalresults["positive for OrgCell"]))))
                  #file="{}/cell_org/hmm-prospector/hmmsearch_results/final.tab".format(param["out"])
                  #positive=get_names(file)
                  if valid and not positive==[]:
                    match=["Invalid","Valid"]
                    if call_fetch(name,o,i,match,positive,True):
                      if os.path.isfile("{}/results/valid_HMMs_list.txt".format(param['out'])):
                        validated.extend(open("{}/results/valid_HMMs_list.txt".format(param['out']),"r").read().split("\n"))
                        for model in validated:
                          if "vHMM_" in model:
                            key="negative for OrgCell"
                            if not key in finalresults.keys():
                              finalresults[key]=[model]
                            elif not name in finalresults[key]:
                              finalresults[key].append(model)
                            log.write("\n{}\n".format(model))
                            log.write("'{}' didn't match cell_org!\n".format(model))
                            log.write("Profile HMM '{}' is valid!\n".format(model))
                        print("{} negative for OrgCell = validated".format(len(set(finalresults["negative for OrgCell"]))))
              elif param['db_type']=='long':
                log.close()
                log=open(os.path.join(param["out"], 'file.log'),'a')
                log.write("\n")
                name='cell_org'
                o=param['out']
                i=param['i']
                match="Match"
                positive=[]
                valid,positive=parse_table(name,o,i,match)
                finalresults["positive for OrgCell"]=positive
                #file="{}/cell_org/hmm-prospector/hmmsearch_results/final.tab".format(param["out"])
                #positive=get_names(file)
                if valid and not positive==[]:
                  match=["Match","Valid"]
                  if call_fetch(name,o,i,match,positive,True):
                    if os.path.isfile("{}/results/valid_HMMs_list.txt".format(param['out'])):
                      validated.extend(open("{}/results/valid_HMMs_list.txt".format(param['out']),"r").read().split("\n"))
                      for model in validated:
                        if "vHMM_" in model:
                          key="negative for OrgCell"
                          if not key in finalresults.keys():
                            finalresults[key]=[model]
                          elif not name in finalresults[key]:
                            finalresults[key].append(model)
                          log.write("\n{}\n".format(model))
                          log.write("'{}' didn't match cell_org!\n".format(model))
                          log.write("Profile HMM '{}' is valid!\n".format(model))
                      print("{} negative for OrgCell = validated".format(len(set(finalresults["negative for OrgCell"]))))
                      print("{} positive for OrgCell".format(len(set(finalresults["positive for OrgCell"]))))
                    log.close()
                    log=open(os.path.join(param["out"], 'file.log'),'a')
                    log.write("\n")
                    table=os.path.join(param['out'],'cell_org','hmm-prospector','table2.csv')
                    folder=param['vir_db']
                    models=os.path.join(param['out'],'cell_org','models')
                    vir_db,error=get_vir_db(table,folder,models)
                    if not vir_db==False:
                      log.close()
                      try:
                        if not os.path.isdir(os.path.join(param['out'],'vir')):
                          os.mkdir(os.path.join(param['out'],'vir'))
                      except Exception as e: 
                        print(e)
                        print("ERROR: vir directory wasn't created!")
                        log.write("\nERROR: vir directory wasn't created!\n")
                      else:
                        log=open(os.path.join(param["out"], 'file.log'),'a')
                        csv=open(os.path.join(param['out'],'scores.csv'),"w")
                        csv.write("name\tV\tC\tP\tC<P\tV-C\tnew_cutoff_score\told_cutoff_score\tnew>old\n")
                        csv.close()
                        for name in sorted(list(vir_db.keys())):
                          csv.close()
                          log.write("\n{}\n".format(name))
                          log.write("'{}' match cell_org!\n".format(name))
                          log.write("'{}' HMM model: {}\n".format(name,os.path.join(param['out'],'cell_org','models','{}.hmm'.format(name))))
                          log.write("'{}' cell_org results: grep '{}' {}\n".format(name,name,os.path.join(param['out'],'cell_org','hmm-prospector','table2.csv')))
                          log.write("'{}' taxon: {}\n".format(name,taxonomy[name]['Taxon']))
                          log.write("'{}' taxon rank: {}\n".format(name,taxonomy[name]['Rank']))
                          if 'Family' in taxonomy[name].keys():
                            log.write("'{}' family: {}\n".format(name,taxonomy[name]['Family']))
                          log.write("'{}' protein: {}\n".format(name,taxonomy[name]['Protein']))
                          log.write("'{}' protein fasta file: {}\n".format(name,vir_db[name]))
                          if not name=="'":
                            try:
                              os.chdir(os.path.join(param['out'],'vir'))
                            except Exception as e: 
                              print(e)
                              print("ERROR: vir directory wasn't opened!")
                              log.write("\nERROR: vir directory wasn't opened!\n")
                            else:
                              db=vir_db[name]
                              i=os.path.join(param['out'],'cell_org','models','{}.hmm'.format(name))
                              o=os.path.join(param['out'],'vir',name)
                              look=os.path.join(param['out'],'vir','{}_prosp_out'.format(name))
                              start=datetime.now()
                              cmd=cmd_prospector(param['reads_vir'],db,i,o,"vir")
                              log.close()
                              log=open(os.path.join(param["out"], 'file.log'),'a')
                              if run_prospector(db,i,o,param['reads_vir'],"vir",look,cmd):
                                #print("'{}' hmm-prospector runtime for {}: {}".format(name,"vir",datetime.now()-start))
                                log.write("'{}' hmm-prospector runtime for {}: {}\n".format(name,"vir",datetime.now()-start))
                                log.write("'{}' vir results: {}\n".format(name,os.path.join(param["out"],"vir",name)))
                                log.close()
                                log=open(os.path.join(param["out"], 'file.log'),'a')
                                try:
                                  if not os.path.isdir(os.path.join(param['out'],'new_scores')):
                                    os.mkdir(os.path.join(param['out'],'new_scores'))
                                  os.chdir(os.path.join(param['out'],'new_scores'))
                                except Exception as e: 
                                  print(e)
                                  print("ERROR: new_scores directory wasn't created!")
                                  log.write("\nERROR: new_scores directory wasn't created!\n")
                                else:
                                  cell_org="{}/cell_org/hmm-prospector/table1.csv".format(param["out"])
                                  vir="{}/vir/{}/table1.csv".format(param["out"],name)
                                  pt=param['pt']
                                  i=os.path.join(param['out'],'cell_org','models','{}.hmm'.format(name))
                                  o=os.path.join(os.path.join(param['out'],'new_scores','{}.hmm'.format(name)))
                                  csv=open(os.path.join(param['out'],'scores.csv'),"a")
                                  if compare_scores(cell_org,vir,pt,i,o,name):
                                    csv.close()
                                    log.close()
                                    log=open(os.path.join(param["out"], 'file.log'),'a')
                                    if len(open(o,"r").read())==0 or not os.path.isfile(o):
                                      log.write("ERROR: Profile HMM file '{}.hmm' with new cutoff score was not created or is empty!\n".format(name))
                                      print("ERROR: Profile HMM file '{}.hmm' with new cutoff score was not created or is empty!".format(name))
                                    else:
                                      log.write("'{}' HMM model with new cutoff score: {}\n".format(name,os.path.join(param['out'],'new_scores','{}.hmm'.format(name))))
                                      try:
                                        if not os.path.isdir(os.path.join(param['out'],'recall')):
                                          os.mkdir(os.path.join(param['out'],'recall'))
                                        os.chdir(os.path.join(param['out'],'recall'))
                                      except Exception as e: 
                                        print(e)
                                        print("ERROR: Recall directory wasn't created!")
                                        log.write("\nERROR: Recall directory wasn't created!\n")
                                      else:
                                        db=vir_db[name]
                                        i=os.path.join(os.path.join(param['out'],'new_scores','{}.hmm'.format(name)))
                                        o=os.path.join(param['out'],'recall',name)
                                        look=os.path.join(param['out'],'recall','{}_prosp_out'.format(name))
                                        start=datetime.now()
                                        cmd=cmd_prospector(param['reads_vir'],db,i,o,"recall")
                                        log.close()
                                        log=open(os.path.join(param["out"], 'file.log'),'a')
                                        if run_prospector(db,i,o,param['reads_vir'],"recall",look,cmd):
                                          #print("'{}' hmm-prospector runtime for {}: {}".format(name,"recall",datetime.now()-start))
                                          log.write("'{}' hmm-prospector runtime for {}: {}\n".format(name,"recall",datetime.now()-start))
                                          log.write("'{}' recall results: {}\n".format(name,os.path.join(param["out"],"recall",name)))
                                          log.close()
                                          log=open(os.path.join(param["out"], 'file.log'),'a')
                                          try:
                                            if not os.path.isdir(os.path.join(param['out'],'results')):
                                              os.mkdir(os.path.join(param['out'],'results'))
                                            os.chdir(os.path.join(param['out'],'results'))
                                          except Exception as e: 
                                            print(e)
                                            print("ERROR: Results directory wasn't created!")
                                            log.write("\nERROR: Results directory wasn't created!\n")
                                          else:
                                            path=param['out']
                                            fasta=vir_db[name]
                                            pd=param['pd']
                                            if validate_recall(path,name,fasta,pd):
                                              validated.append(name)
                                              try:
                                                if not os.path.isdir(os.path.join(param['out'],'results','valid')):
                                                  os.mkdir(os.path.join(param['out'],'results','valid'))
                                                os.chdir(os.path.join(param['out'],'results','valid'))
                                              except Exception as e: 
                                                print(e)
                                                print("ERROR: Valid directory wasn't created!")
                                                log.write("\nERROR: Valid directory wasn't created!\n")
                                              else:
                                                try:
                                                  os.system('cp {}/new_scores/{}.hmm {}/results/valid/{}.hmm'.format(param['out'],name,param['out'],name))
                                                except Exception as e: 
                                                  print(e)
                                                  print("ERROR: {} file cannot be copied to valid folder!".format(name))
                                                  log.write("ERROR: {} file cannot be copied to valid folder!\n".format(name))
                                                else:
                                                  #print("Profile HMM '{}' is valid!".format(name))
                                                  log.write("Profile HMM '{}' is valid!\n".format(name))
                        if "invalidated because the highest Cscore > lowest Vscore" in finalresults.keys():
                          print("{} invalidated because the highest Cscore > lowest Vscore".format(len(set(finalresults["invalidated because the highest Cscore > lowest Vscore"]))))
                        if "invalidated because the highest Cscore > (pt%) x lowest Vscore" in finalresults.keys():
                          print("{} invalidated because the highest Cscore > (pt%) x lowest Vscore".format(len(set(finalresults["invalidated because the highest Cscore > (pt%) x lowest Vscore"]))))
                        if "invalidated because recall < (pd%)" in finalresults.keys():
                          print("{} invalidated because recall < (pd%)".format(len(set(finalresults["invalidated because recall < (pd%)"]))))
                        if "validated with adjusted cutoff score" in finalresults.keys():
                          print("{} validated with adjusted cutoff score".format(len(set(finalresults["validated with adjusted cutoff score"]))))
                log.write("\n")
                #print("\n")
                if not os.path.isdir(os.path.join(param['out'],'results')):
                  os.mkdir(os.path.join(param['out'],'results'))
                os.chdir(os.path.join(param['out'],'results'))
                if validated==[]: 
                  print("ERROR: List of valid HMMs is empty!")
                  log.write("ERROR: List of valid HMMs is empty!\n")
                else:
                  log.write("\n{} valid HMMs: {}\n".format(len(validated),validated))
                  #vfile=open("{}/results/validHMMs".format(param['out']),"w")
                  #vfile.write("\n".join(validated))
                  #vfile.close()
                  if os.path.isfile("{}/results/valid/valid.hmm".format(param['out'])):
                    os.remove("{}/results/valid/valid.hmm".format(param['out']))
                  try:
                    subprocess.check_output("cat {}/results/valid/*.hmm > {}/results/valid/valid.hmm".format(param['out'],param['out']), shell=True)
                  except:
                    print("ERROR: All {} valid HMM were not concatenated!".format(len(validated)))
                    log.write("ERROR: All {} valid HMM were not concatenated!\n".format(len(validated)))
                  else:
                    print("All {} valid HMM were concatenated!".format(len(validated)))
                    log.write("All {} valid HMM were concatenated!\n".format(len(validated)))
                i=param['i']
                negative=[]
                all=validated
                all.extend(error)
                valid,negative=get_diff(i,all,"Valid","")
                name=""
                o=param['out']
                match=["Invalid"]
                both=False   
                if negative==[]:
                  print("ERROR: List of invalid HMMs is empty!")
                  log.write("ERROR: List of invalid HMMs is empty!\n")
                else:
                  #log.write("\nInvalid: {}\n".format(negative))     
                  call_fetch(name,o,i,match,negative,both)
                error=sorted(list(set(error)))
                log.write("\n{} Error HMMs: {}\n".format(len(error),error))
                total=int(str(subprocess.check_output("grep 'NAME' {}| wc -l".format(param['i']), shell=True)).replace("b'","").replace("\\n'","").split()[0])
                if len(positive)+len(negative)+len(error)==total:
                  print("ERROR: Something is wrong, the sum of invalid and valid models is not equal to the total number of models!")
                  log.write("ERROR: Something is wrong, the sum of invalid and valid models is not equal to the total number of models!\n")
              log.write("\nFinal result:\n")
              total=int(str(subprocess.check_output("grep 'NAME' {}| wc -l".format(param['i']), shell=True)).replace("b'","").replace("\\n'","").split()[0])
              log.write("{} tested profile HMMs\n".format(total))
              if param['db_type']=='short': 
                log.write("\t{} positive for OrgCell = invalidated\n".format(len(set(finalresults["positive for OrgCell"]))))
              elif param['db_type']=='long':
                log.write("\t{} positive for OrgCell\n".format(len(set(finalresults["positive for OrgCell"]))))
                if "invalidated because the highest Cscore > lowest Vscore" in finalresults.keys():
                  log.write("\t\t{} invalidated because the highest Cscore > lowest Vscore\n".format(len(set(finalresults["invalidated because the highest Cscore > lowest Vscore"]))))
                if "invalidated because the highest Cscore > (pt%) x lowest Vscore" in finalresults.keys():
                  log.write("\t\t{} invalidated because the highest Cscore > (pt%) x lowest Vscore\n".format(len(set(finalresults["invalidated because the highest Cscore > (pt%) x lowest Vscore"]))))
                if "invalidated because recall < (pd%)" in finalresults.keys():
                  log.write("\t\t{} invalidated because recall < (pd%)\n".format(len(set(finalresults["invalidated because recall < (pd%)"]))))
                if "validated with adjusted cutoff score" in finalresults.keys():
                  log.write("\t\t{} validated with adjusted cutoff score\n".format(len(set(finalresults["validated with adjusted cutoff score"]))))
              log.write("\t{} negative for OrgCell = validated\n".format(len(set(finalresults["negative for OrgCell"]))))
            execution=datetime.now() - start_time
            print("\nExecution time: {}".format(execution))
            log.write("\nExecution time: {}".format(execution))
            log.close()
print("Done.")
quit()