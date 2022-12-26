#!/usr/bin/env python3

version='1.0.1'

import os
from datetime import datetime
import argparse
from argparse import RawTextHelpFormatter
import sys
import subprocess

param=dict()
validated=[]
invalidated=[]
start_time = datetime.now()
call=os.path.abspath(os.getcwd())

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
  
  string = '(c) 2022. Arthur Gruber & Giuliana Pola\n'
  string += 'For more information access: https://github.com/GiulianaPola/model_validator\n'
  string += "\nUsage:\n"
  string += "model_validator.py -conf <configuration_file>\n"
  string += "model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'short' -cell_org_db_frag <fasta_file> <optional parameters>\n"
  string += "model_validator.py -i <hmm_file> -model_type <'l'|'long'|'s'|'short'> -db_type 'long' -cell_org_db <fasta_file> -vir_db <directory> <optional parameters>\n"
  string += "\nMandatory parameters:\n"
  
  size = int(str(os.get_terminal_size()).split("columns=")[-1].split(",")[0])
  #print(size)
  
  keys=list(mandatory.keys())
  keys.extend(list(optional.keys()))
  values=list(mandatory.values())
  values.extend(list(optional.values()))
  maxk=0
  maxv=0
  
  for key in keys:
    if len(key)>maxk:
      maxk=len(key)
  
  rest=size-maxk-1
  
  for key in sorted(mandatory.keys()):
    line=''
    line+='{message: <{width}} '.format(
    message=key,
    fill=' ',
    align='<',
    width=maxk,
    )
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
      newname=os.path.join(path, str(name+str(i)))
  elif typ=='file':
    while os.path.isfile(newname):
      i+=1
      newname=os.path.join(path, str(name+str(i)))
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
              for line in text.split('\n'):
                  index=line.find("=")
                  value=line[index+1:]
                  key=line[0:index].lower()
                  if key=='input_file' or key=='i':
                      args['i']=value.replace('"','').strip()
                  elif key=='output' or key=='out' or key=='o':
                      args['out']=value
                  elif key=='db_type' or key=='dbtype':
                      args['db_type']=value
                  elif key=='modeltype' or key=='model_type':
                      args['model_type']=value
                  elif key=='cell_org_db':
                      args['cell_org_db']=value.replace('"','').strip()
                  elif key=='cell_org_db_frag':
                      args['cell_org_db_frag']=value.replace('"','').strip()
                  elif key=='vir_db':
                      args['vir_db']=value.replace('"','').strip()
                  elif key=='vir_db_frag':
                      args['vir_db_frag']=value
                  elif key=='pt':
                      args['pt']=value
                  elif key=='pd':
                      args['pd']=value
                  elif key in ['long_hmms_x_short_reads_cell','lh_sr']:
                      args['long_hmms_x_short_reads_cell']=value
                  elif key in ['short_hmms_x_short_reads_cell','sh_sr']:
                      args['short_hmms_x_short_reads_cell']=value
                  elif key in ['long_hmms_x_long_reads_cell','lh_lr']:
                      args['long_hmms_x_long_reads_cell']=value
                  elif key in ['short_hmms_x_long_reads_cell','sh_lr']:
                      args['short_hmms_x_long_reads_cell']=value
                  elif key in ['long_hmms_x_reads_vir','lh_rv']:
                      args['long_hmms_x_reads_vir']=value
                  elif key in ['short_hmms_x_reads_vir','sh_rv']:
                      args['short_hmms_x_reads_vir']=value
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
              args['short_hmms_x_long_reads_cell']="'-cpu 20 -n 10 -r no -sv 1.0'"
              args['short_hmms_x_reads_vir']="'-cpu 20 -n 10 -r no -sv 1.0'"
              param['reads_cell']="-cpu 20 -n 10 -r no -sv 1.0" 
              param['reads_vir']="-cpu 20 -n 10 -r no -sv 1.0"
            elif param['model_type']=='long':
              args['long_hmms_x_long_reads_cell']="'-cpu 20 -n 10 -r no -sv 0.6'" 
              args['long_hmms_x_reads_vir']="'-cpu 20 -n 10 -r no -sv 0.6'" 
              param['reads_cell']="-cpu 20 -n 10 -r no -sv 0.6"
              param['reads_vir']="-cpu 20 -n 10 -r no -sv 0.6"
        elif param['db_type']=='short':
          if 'reads_cell' not in param.keys():
            if param['model_type']=='short':
              args['short_hmms_x_short_reads_cell']="'-cpu 20 -n 10 -r no -rl c -sv 1.5'"
              param['reads_cell']="-cpu 20 -n 10 -r no -rl c -sv 1.5"
            elif param['model_type']=='long':
              args['long_hmms_x_short_reads_cell']="'-cpu 20 -n 10 -r no -rl c -sv 0.8'"
              param['reads_cell']="-cpu 20 -n 10 -r no -rl c -sv 0.8"

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

def run_prospector(db,i,o,parameters,name,look):
  positive=[]
  valid=True
  cmd="nohup hmm-prospector.pl"
  cmd+=" -d {}".format(db)
  cmd+=" -i {}".format(i)
  if "hmm-prospector.pl" in parameters:
    parameters=parameters.split('hmm-prospector.pl')[1]
  for arg in parameters.split(' -'):
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
  log.write("{}\n".format(cmd))
  start=datetime.now()
  os.system(cmd)
  while not os.path.isfile(look):
    pass
  output=open(look,"r")
  process=os.popen("ps -u {} --sort +start_time".format(os.getlogin())).read()
  table="{}/hmm-prospector/table2.csv".format(o)
  while "hmm-prospector." in process and "Done." not in output.read() and not os.path.isfile(table):
    output.close()
    output=open(look,"r")
    process=os.popen("ps -u {} --sort +start_time".format(os.getlogin())).read()
  output.close()
  log.write("hmm-prospector runtime for {}: {}\n".format(name,datetime.now()-start))
  if name=="cell_org":
    print("hmm-prospector runtime for {}: {}".format(name,datetime.now()-start))
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
    mismatch="Mismatch"
  elif match=="Mismatch":
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
        listHMMs=os.path.join(o,name,"matchedHMMs")
      elif match=="Valid" or match=="Invalid":
        if not os.path.isdir(os.path.join(o,"results")):
          os.mkdir(os.path.join(o,"results"))
        listHMMs=os.path.join(o,"results","{}HMMs".format(match.lower()))
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
            selectedfile=open(os.path.join(o, 'selected.hmm'),'w')
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
    mismatch="Mismatch"
  elif match[0]=="Mismatch":
    mismatch="Match"
  if not name=="":
    mismatch+=" for {}".format(name)
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
    if not fetch_hmm(positive,i,out):
      valid=False
      if match[0] in ["Valid","Invalid"]:
        print("{} HMMs couldn't be fetched!".format(match[0].capitalize()))
        log.write("{} HMMs couldn't be fetched!\n".format(match[0].capitalize()))
      else:
        print("HMMs matching {} couldn't be fetched!".format(name))
        log.write("HMMs matching {} couldn't be fetched!\n".format(name))
    else:
      if match[0] =="Match":
        log.write("Runtime to fetch the HMMs matching {}: {}\n".format(name,datetime.now()-start))
        print("{} HMMs matching {} were fetched!".format(len(positive),name))
        log.write("{} HMMs matching {} were fetched!\n".format(len(positive),name))
        log.write("{} match HMMs: {}\n".format(name,positive))
      elif match[0]=="Invalid" or match[0]=="Valid":
        log.write("{} HMMs: {}\n".format(match[0].capitalize(),positive))
        print("Runtime to fetch {} HMMs: {}".format(match[0].lower(),datetime.now()-start))
        log.write("Runtime to fetch {} HMMs: {}\n".format(match[0].lower(),datetime.now()-start))
        print("{} {} HMMs were fetched!".format(len(positive),match[0].lower()))
        log.write("{} {} HMMs were fetched!\n".format(len(positive),match[0].lower()))
  if both==True:
    valid,negative=get_diff(i,positive,match[0],name)
    if negative==[]:
      log.write("{} HMMs list is empty!\n".format(mismatch))
    else:
        if not os.path.isdir(os.path.join(o,"results")):
          try:
            os.mkdir(os.path.join(o,"results"))
          except:
            print("'results' directory was not created!")
            log.write("'results' directory was not created!\n")
            valid=False
          else:
            print("'results' directory was created!")
            log.write("'results' directory was created!\n")
        log.write("{}:{}\n".format(mismatch,negative))
        valid,listHMMs=list_HMMs("",match[1],o,negative)
        if valid and not listHMMs=='':
          if match[1]=="Valid" or match[1]=="Invalid":
            out=os.path.join(o,"results",match[1].lower())
            start=datetime.now()
            if not fetch_hmm(negative,i,out):
              valid=False
              if not match[1]==mismatch.split(" ")[0]:
                print("{} HMMs couldn't be fetched!".format(mismatch.capitalize()))
                log.write("{} HMMs couldn't br fetched!\n".format(mismatch.capitalize()))
              else:
                print("{} HMMs couldn't be fetched!".format(match[1].capitalize()))
                log.write("{} HMMs couldn't be fetched!\n".format(match[1].capitalize()))
            else:
              if not match[1]==mismatch.split(" ")[0]:
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

def get_rank(txid,taxon):
  import requests
  link = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}&lvl=3&lin=f&keep=1&srchmode=1&unlock".format(txid)
  try:
    html = requests.get(link).text
  except Exception as e: 
    log.write("{}\n".format(e))
    print("Rank for taxon {}({}) was not found!".format(taxon,txid))
    log.write("Rank for taxon {}({}) was not found!\n".format(taxon,txid))
    valid=False
  else:
    if "family" in html:
      valid="family"
    elif "genus" in html:
      valid="genus"
    else:
      valid=False
  if valid==False:
    print("Rank for taxon {}({}) was not found!".format(taxon,txid))
    log.write("Rank for taxon {}({}) was not found!\n".format(taxon,txid))
  return valid

def get_vir_db(models,path):
  search=dict()
  taxons=[]
  vir_db=dict()
  proteins=dict()
  tax_dir=''
  prot_dir=''
  taxonomy=dict()
  for model in open(models,"r").read().split("//"):
    if 'NAME' in model and 'TAXON' in model and 'TXID' in model and 'PROTEIN' in model:
      name=model.split("NAME")[-1].split("\n")[0].strip()
      taxon=model.split("TAXON")[-1].split("\n")[0].strip()
      if taxon not in taxons and not taxon=="'":
        taxons.append(taxon)
      txid=model.split("TXID")[-1].split("\n")[0].strip()
      if taxon not in taxonomy.keys() and not taxon=="'":
        taxonomy[taxon]=txid
      protein=model.split("PROTEIN")[-1].split("\n")[0].strip()
      if not taxon in proteins.keys() and not taxon=="'":
        proteins[taxon]=[protein]
      else:
        proteins[taxon].append(protein)
      if "{}/{}".format(taxon,protein) in search.keys() and not "'" in [taxon,protein,name]:
        search["{}/{}".format(taxon,protein)].append(name)
      else:
        search["{}/{}".format(taxon,protein)]=[name]
  for taxon in set(taxons):
    rank=get_rank(taxonomy[taxon],taxon)
    if not rank==False:
      tax_dir=get_path(path,taxon)
      if tax_dir==False:
        print("Missing taxon '{}' folder!".format(taxon))
        log.write("Missing taxon '{}' folder!\n".format(taxon))
      else:
        for protein in set(proteins[taxon]):
          prot_file=''
          prot_dir=get_path(path,protein)
          if prot_dir==False:
            n=0
            for filename in os.listdir(tax_dir):
              if os.path.isdir(os.path.join(tax_dir,filename)) and filename.replace(" ","_").lower() in protein.replace(" ","_").lower():
                if len(filename)>n:
                  prot_dir=os.path.join(tax_dir,filename)
                  n=len(filename)
              if os.path.isdir(os.path.join(tax_dir,filename)) and protein.replace(" ","_").lower() in filename.replace(" ","_").lower():
                if len(protein)>n:
                  prot_dir=os.path.join(tax_dir,filename)
                  n=len(protein)
          if prot_dir==False:
            print("Missing protein '{}' folder in '{}' folder!".format(protein,taxon))
            log.write("Missing protein '{}' folder in '{}' folder!\n".format(protein,taxon))
          else:
            if not os.path.isdir(os.path.join(prot_dir,"fasta_sequences")):
              print("Missing 'fasta_sequences' folder in {}/{}!".format(protein,taxon))
              log.write("Missing 'fasta_sequences' folder in {}/{}!\n".format(protein,taxon))
            else:
              for name in search["{}/{}".format(taxon,protein)]:
                if not name=="'":
                  if rank=='family':
                    prot_file=get_path(os.path.join(prot_dir,"fasta_sequences"),"{}_no_redundancy.fasta".format(taxon))
                    if prot_file==False:
                      print("Missing protein '{}_no_redundancy.fasta' file in {}/{}/fasta_sequences!".format(taxon,taxon,protein))
                      log.write("Missing protein '{}_no_redundancy.fasta' file in {}/{}/fasta_sequences!\n".format(taxon,taxon,protein))
                    else:
                      vir_db[name]=prot_file
                  if rank=='genus':
                    prot_file=get_path(os.path.join(prot_dir,"fasta_sequences"),"{}.fasta".format(taxon))
                    if prot_file==False:
                      print("Missing protein '{}.fasta' file in {}/{}/fasta_sequences!".format(taxon,taxon,protein))
                      log.write("Missing protein '{}.fasta' file in {}/{}/fasta_sequences!\n".format(taxon,taxon,protein))
                    else:
                      vir_db[name]=prot_file
  if len(vir_db)==0:
    return False
  else:
    return vir_db

def get_score(name,file,mode):
  scores=[]
  lines=str(subprocess.check_output("grep ' {} ' {}".format(name,file), shell=True)).split("\\n")
  if "'" in lines:
    lines.remove("'")
  for line in lines:
    columns=line.split(" ")
    while '' in columns:
      columns.remove('')
    while '-' in columns:
      columns.remove('-')
    if len(columns)>=4:
      scores.append(columns[3])
  if mode=="min":
    return min(scores)
  elif mode=="max":
    return min(scores)

def compare_scores(cell_org,vir,pt,i,o,name):
  valid=True
  import subprocess
  headers = str(subprocess.check_output("head {} -n 1".format(vir), shell=True)).replace("b'","").replace("\\n","").split("\\t")
  indexes=[]
  indexes.append(headers.index("query_pHMM")+1)
  indexes.append(headers.index("Score")+1)
  vir_hmm=open(i,"r").read()
  new_hmm=''
#  file="{}/vir/{}/hmmsearch_results/final.tab".format(param["out"],name)
#  mode="min"
#  score=get_score(name,file,mode)
  try:
    V = float(str(subprocess.check_output("grep '{}' {} | cut -f{} | sort | head -n 1".format(name,vir,indexes[1]), shell=True)).replace("b'","").replace("\\n'",""))  
    #V = float(score)
  except Exception as e: 
    log.write("All '{}' vir matches were bellow cutoff score!\n".format(name))
    print("All '{}' vir matches were bellow cutoff score!".format(name))
    valid=False
  else:
    #file="{}/cell_org/hmm-prospector/hmmsearch_results/final.tab".format(param["out"])
    #mode="max"
    #score=get_score(name,file,mode)
    try:
      C = float(str(subprocess.check_output("grep '{}' {} | cut -f{} | sort -rn | head -n 1".format(name,cell_org,indexes[1]), shell=True)).replace("b'","").replace("\\n'",""))
      #C = float(score)
    except Exception as e: 
      log.write("All '{}' cell_org matches were bellow cutoff score!\n".format(name))
      print("All '{}' cell_org matches were bellow cutoff score!".format(name))
      valid=False
    else:
      P = float(V) * pt/100
      if not C < P:
        log.write("'{}' maximum score for cell_org is less than or equal to his minimum score for vir!\n".format(name))
        print("'{}' maximum score for cell_org is less than or equal to his minimum score for vir!".format(name))
        log.write("query_pHMM={}\tmin_vir_score={:.1f}\tmax_cell_org_score={:.1f}\tpt%_min_vir_score={:.1f}\n".format(name,float(V),float(C),float(P)))
        valid=False
      else:
        cutoff = (pt/100 * (V - C)) + C
        log.write("query_pHMM={}\tmin_vir_score={:.1f}\tmax_cell_org_score={:.1f}\t_min_vir_score={:.1f}\tnew_cutoff_score={:.1f}\n".format(name,float(V),float(C),float(P),cutoff))
        if vir_hmm=='':
          log.write("Profile HMM file '{}.hmm' with new cutoff score was not created!\n".format(name))
          print("Profile HMM file '{}.hmm' with new cutoff score was not created!".format(name))
          valid=False
        else:
          new_hmm+=vir_hmm.split("CUTOFF SCORE")[0]+"CUTOFF SCORE    "+str(cutoff)+"\nTAXON"+vir_hmm.split("TAXON")[-1]
          try:
            out=open(o,"w")
            out.write("")
          except:
            log.write("Profile HMM file '{}.hmm' with new cutoff score was not created!\n".format(name))
            print("Profile HMM file '{}.hmm' with new cutoff score was not created!".format(name))
            valid=False
          else:
            out.write(new_hmm)
            out.close()
  return valid
  
def validate_recall(path,name,fasta,pd):
  valid=True
  match=int(str(subprocess.check_output("wc -l {}/recall/{}/table1.csv".format(path,name), shell=True)).replace("b'","").split()[0])-1
  if match<=0:
    print("'{}' recall didn't have any matches!")
    log.write("'{}' recall didn't have any matches!\n")
    valid=False
  else:
    total=int(str(subprocess.check_output("grep '>' {}| wc -l".format(fasta), shell=True)).replace("b'","").replace("\\n'","").split()[0])
    if match/total*100<pd:
      print("'{}' recall detection percentage {:.1f}% is less than {:.1f}% (-pd)!".format(name,match/total*100,pd))
      log.write("'{}' recall detection percentage {:.1f}% is less than {:.1f}% (-pd)!\n".format(name,match/total*100,pd))
      valid=False
    else:
      log.write("'{}' recall detection percentage {:.1f}% is higher than {:.1f}% (-pd)!\n".format(name,match/total*100,pd))
  return valid 

if not len(sys.argv)>1:
    print_help()
elif arguments.h == True:
    print_help()
elif arguments.version == True:
    print(version) 
else:
    if arguments.conf==None:
        args={k: v for k, v in vars(arguments).items() if v is not None and v is not True and v is not False}     
    else:
        valid,args=validate_conf(arguments.conf)
    if valid:
      valid,param=validate_args(args)
      if valid==False:
        print("Input parameters are invalid!")
      if valid==True:
        print("Input parameters are valid!")
        os.chdir(param['out'])
        try:
          log=open(os.path.join(param["out"], 'file.log'),'w')
          log.write('model_validator v{}\n'.format(version))
        except Exception as e: 
          print(e)
          print('Log file was not created!')
        else:
          log.write('\nDatetime:\n{}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
          
          log.write('\nWorking directory:\n{}\n'.format(call))
          
          log.write('\nCommand line:\n{}\n'.format(' '.join(sys.argv)))
          
          log.write('\nUser:\n{}\n'.format(os.getlogin()))
          
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
            print("cell_org directory wasn't created!")
            log.write("\ncell_org directory wasn't created!\n")
          else:
            db=param['cell_org_db']
            i=param['i']
            o=os.path.join(param['out'],"cell_org")
            parameters=param['reads_cell']
            name='cell_org'
            look=os.path.join(param['out'],'cell_org','cell_org_prosp_out')
            if run_prospector(db,i,o,parameters,name,look):
              log.close()
              log=open(os.path.join(param["out"], 'file.log'),'a')
              if param['db_type']=='short':
                try:
                  if not os.path.isdir(os.path.join(param['out'],'results')):
                    os.mkdir(os.path.join(param['out'],'results'))
                  os.chdir(os.path.join(param['out'],'results'))
                except Exception as e: 
                  log.write(e)
                  print("\nResults directory wasn't created!")
                  log.write("\nResults directory wasn't created!\n")
                else:
                  print("Results directory was created!")
                  log.write("\nResults directory was created!\n")
                  name='cell_org'
                  o=param['out']
                  i=param['i']
                  match="Invalid"
                  positive=[]
                  valid,positive=parse_table(name,o,i,match)
                  #file="{}/cell_org/hmm-prospector/hmmsearch_results/final.tab".format(param["out"])
                  #positive=get_names(file)
                  if valid and not positive==[]:
                    match=["Invalid","Valid"]
                    call_fetch(name,o,i,match,positive,True)
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
                #file="{}/cell_org/hmm-prospector/hmmsearch_results/final.tab".format(param["out"])
                #positive=get_names(file)
                if valid and not positive==[]:
                  match=["Match","Valid"]
                  if call_fetch(name,o,i,match,positive,True):
                    if os.path.isfile("{}/results/validHMMs".format(param['out'])):
                      validated.extend(open("{}/results/validHMMs".format(param['out']),"r").read().split("\n"))
                    log.close()
                    log=open(os.path.join(param["out"], 'file.log'),'a')
                    log.write("\n")
                    models=os.path.join(param['out'],'cell_org','models','selected.hmm')
                    path=param['vir_db']
                    vir_db=get_vir_db(models,path)
                    if not vir_db==False:
                      log.close()
                      try:
                        if not os.path.isdir(os.path.join(param['out'],'vir')):
                          os.mkdir(os.path.join(param['out'],'vir'))
                      except Exception as e: 
                        print(e)
                        print("vir directory wasn't created!")
                        log.write("\nvir directory wasn't created!\n")
                      else:
                        log=open(os.path.join(param["out"], 'file.log'),'a')
                        print("\n")
                        for name in sorted(list(vir_db.keys())):
                          log.write("\n{}\n".format(name))
                          if not name=="'":
                            try:
                              os.chdir(os.path.join(param['out'],'vir'))
                            except Exception as e: 
                              print(e)
                              print("vir directory wasn't opened!")
                              log.write("\nvir directory wasn't opened!\n")
                            else:
                              db=vir_db[name]
                              i=os.path.join(param['out'],'cell_org','models','{}.hmm'.format(name))
                              o=os.path.join(param['out'],'vir',name)
                              parameters=param['reads_vir']
                              look=os.path.join(param['out'],'vir','{}_prosp_out'.format(name))
                              if run_prospector(db,i,o,parameters,"vir",look):
                                log.close()
                                log=open(os.path.join(param["out"], 'file.log'),'a')
                                try:
                                  if not os.path.isdir(os.path.join(param['out'],'new_scores')):
                                    os.mkdir(os.path.join(param['out'],'new_scores'))
                                  os.chdir(os.path.join(param['out'],'new_scores'))
                                except Exception as e: 
                                  print(e)
                                  print("new_scores directory wasn't created!")
                                  log.write("\nnew_scores directory wasn't created!\n")
                                else:
                                  cell_org="{}/cell_org/hmm-prospector/table1.csv".format(param["out"])
                                  vir="{}/vir/{}/table1.csv".format(param["out"],name)
                                  pt=param['pt']
                                  i=os.path.join(param['out'],'cell_org','models','{}.hmm'.format(name))
                                  o=os.path.join(os.path.join(param['out'],'new_scores','{}.hmm'.format(name)))
                                  if compare_scores(cell_org,vir,pt,i,o,name):
                                    log.close()
                                    log=open(os.path.join(param["out"], 'file.log'),'a')
                                    if len(open(o,"r").read())==0 or not os.path.isfile(o):
                                      log.write("Profile HMM file '{}.hmm' with new cutoff score was not created or is empty!\n".format(name))
                                      print("Profile HMM file '{}.hmm' with new cutoff score was not created or is empty!".format(name))
                                    else:
                                      try:
                                        if not os.path.isdir(os.path.join(param['out'],'recall')):
                                          os.mkdir(os.path.join(param['out'],'recall'))
                                        os.chdir(os.path.join(param['out'],'recall'))
                                      except Exception as e: 
                                        print(e)
                                        print("Recall directory wasn't created!")
                                        log.write("\nRecall directory wasn't created!\n")
                                      else:
                                        db=vir_db[name]
                                        i=os.path.join(os.path.join(param['out'],'new_scores','{}.hmm'.format(name)))
                                        o=os.path.join(param['out'],'recall',name)
                                        parameters=param['reads_vir']
                                        look=os.path.join(param['out'],'recall','{}_prosp_out'.format(name))
                                        if run_prospector(db,i,o,parameters,"recall",look):
                                          log.close()
                                          log=open(os.path.join(param["out"], 'file.log'),'a')
                                          try:
                                            if not os.path.isdir(os.path.join(param['out'],'results')):
                                              os.mkdir(os.path.join(param['out'],'results'))
                                            os.chdir(os.path.join(param['out'],'results'))
                                          except Exception as e: 
                                            print(e)
                                            print("Results directory wasn't created!")
                                            log.write("\nResults directory wasn't created!\n")
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
                                                print("Valid directory wasn't created!")
                                                log.write("\nValid directory wasn't created!\n")
                                              else:
                                                try:
                                                  os.system('cp {}/new_scores/{}.hmm {}/results/valid/{}.hmm'.format(param['out'],name,param['out'],name))
                                                except Exception as e: 
                                                  print(e)
                                                  print("{} file cannot be copied to valid folder!".format(name))
                                                  log.write("{} file cannot be copied to valid folder!\n".format(name))
                                                else:
                                                  print("Profile HMM '{}' is valid!".format(name))
                                                  log.write("Profile HMM '{}' is valid!\n".format(name))
                log.write("\n")
                print("\n")
                if not os.path.isdir(os.path.join(param['out'],'results')):
                  os.mkdir(os.path.join(param['out'],'results'))
                os.chdir(os.path.join(param['out'],'results'))
                if validated==[]: 
                  print("List of valid HMMs is empty!")
                  log.write("List of valid HMMs is empty!\n")
                else:
                  log.write("\nValid: {}\n".format(len(validated),validated))
                  vfile=open("{}/results/validHMMs".format(param['out']),"w")
                  vfile.write("\n".join(validated))
                  vfile.close()
                  if os.path.isfile("{}/results/valid/selected.hmm".format(param['out'])):
                    os.remove("{}/results/valid/selected.hmm".format(param['out']))
                  try:
                    subprocess.check_output("cat {}/results/valid/*.hmm > {}/results/valid/selected.hmm".format(param['out'],param['out']), shell=True)
                  except:
                    print("All {} valid HMM were not concatenated!".format(len(validated)))
                    log.write("All {} valid HMM were not concatenated!\n".format(len(validated)))
                  else:
                    print("All {} valid HMM were concatenated!".format(len(validated)))
                    log.write("All {} valid HMM were concatenated!\n".format(len(validated)))
                i=param['i']
                negative=[]
                valid,negative=get_diff(i,validated,"Valid","")
                name=""
                o=param['out']
                match=["Invalid"]
                both=False   
                if negative==[]:
                  print("List of invalid HMMs is empty!")
                  log.write("List of invalid HMMs is empty!\n")
                else:
                  #log.write("\nInvalid: {}\n".format(negative))     
                  call_fetch(name,o,i,match,negative,both)
                total=int(str(subprocess.check_output("grep 'NAME' {}| wc -l".format(param['i']), shell=True)).replace("b'","").replace("\\n'","").split()[0])
                if len(positive)+len(negative)==total:
                  print("Something is wrong, the sum of invalid and valid models is not equal to the total number of models!")
                  log.write("Something is wrong, the sum of invalid and valid models is not equal to the total number of models!\n")
            execution=datetime.now() - start_time
            print("\nExecution time: {}".format(execution))
            log.write("\nExecution time: {}".format(execution))
            log.close()
print("Done.")