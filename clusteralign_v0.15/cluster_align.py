# takes input files, splits into pieces, copy to all nodes and run analysis
import sys
import os
import shutil
import time as ztime
import subprocess 
import string
from random import choice
import datetime
import pprint
from datetime import datetime, date, time

"""
clusteralign.py

Created by Warren Emmett on 2010-03-22.
Copyright (c) 2010 Warren Emmett. All rights reserved.
"""


########################
class progressBar:
	def __init__(self, minValue = 0, maxValue = 10, totalWidth=20):
		self.progBar = "[]"   # This holds the progress bar string
		self.min = minValue
		self.max = maxValue
		self.span = maxValue - minValue
		self.width = totalWidth
		self.amount = 0       # When amount == max, we are 100% done 
		self.updateAmount(0)  # Build progress bar string

	def updateAmount(self, newAmount = 0):
		if newAmount < self.min: newAmount = self.min
		if newAmount > self.max: newAmount = self.max
		self.amount = newAmount

		# Figure out the new percent done, round to an integer
		diffFromMin = float(self.amount - self.min)
		percentDone = (diffFromMin / float(self.span)) * 100.0
		percentDone = round(percentDone)
		percentDone = int(percentDone)

		# Figure out how many hash bars the percentage should be
		allFull = self.width - 2
		numHashes = (percentDone / 100.0) * allFull
		numHashes = int(round(numHashes))

		# build a progress bar with hashes and spaces
		self.progBar = "[" + '#'*numHashes + ' '*(allFull-numHashes) + "]"

		# figure out where to put the percentage, roughly centered
		percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
		percentString = str(percentDone) + "%"

		# slice the percentage into the bar
		self.progBar = self.progBar[0:percentPlace] + percentString + self.progBar[percentPlace+len(percentString):]

	def __str__(self):
		return str(self.progBar)

#########################
def collect_stats(job_list,cluster_data,run_data,tag):
    # collect stats on each job and write to file
	stats_list = []
	
	jobs_done = []
	
	while 1:
	    
		for job in job_list:
			
					
			if job in jobs_done:
				#print job
				#print jobs_done
				#print job_list
				#raw_input()
				continue
			
			else:
				#print "doing",job
				#raw_input()
				stats_proc = subprocess.Popen([cluster_data["cluster_acct"],cluster_data["cluster_qacct_us"],run_data["username"],cluster_data["cluster_js"],job],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			    
				for line in stats_proc.stdout:
				
					if line.startswith("error:"):
						
						break
					else:
						stats_list.append(line)
				
			    
				jobs_done.append(job)   
				
		if len(jobs_done) == len(job_list):
				break
	    
			
		    
			 
	open(os.path.join(aligner_data[aligner]["output_directory"],"%s_stats_%s_run.txt"%(tag,aligner)),"w").writelines(stats_list)

def delete_tmp_dirs(run_dirs,run_data,diag):
	# delete all temporary files etc
	#shutil.rmtree(run_dirs["script_dir"])
	shutil.rmtree(run_dirs["dir_for_read_parts"])
	shutil.rmtree(run_dirs["tmp_results_dir"])
	shutil.rmtree(run_dirs["dir_job_done"])
	#if diag == 1:
		

######################################################################################
def split_file(file_to_split,desired_file_size,dest,not_strict):

    
	file_root = os.path.split(file_to_split)[-1].split(".")[0]
	
	
	
	total_lines = 0
	for line in open(file_to_split):
		total_lines+=1
	
	#"""
	if not_strict == False:
		if total_lines%4 > 0:
			print >> sys.stderr,"Inconsistent filesize! Please remove comments from file"
			sys.exit()
	#"""         
	file_size= os.path.getsize(file_to_split)/1024/1024
	       
	split_number = int(4*desired_file_size)    
	print >> sys.stderr,"Number of lines per file: %s"%(split_number)
	
	log_list.append(str(file_size)+"\t"+str(desired_file_size)+"\t"+str(split_number)+"\n")
	
	if split_number < total_lines:
		file_list = []
		outlist = []    
		count = 1
		tot_c = 0
		for line in open(file_to_split):
			tot_c +=1
			outlist.append(line)
			
			if tot_c >= split_number:
			    
				if line[0] in ["A","C","T","G","N","@","+"]:
					continue
				file_name = os.path.join(dest,file_root+"PART_"+str(count+10))
				file_list.append(file_name)
				open(os.path.join(dest,file_name),"w").writelines(outlist)
				outlist = []
				count +=1
				tot_c = 0
		
		if len(outlist) > 0:
			file_list.append(os.path.join(dest,file_root+"PART_"+str(count+10)))
			open(os.path.join(dest,file_root+"PART_"+str(count+10)),"w").writelines(outlist)
			outlist = []  
		
		return file_list
	else:
		return [file_to_split]
	
	

def paired_split_files(read_file_name1,read_file_name2,divide_by,run_dirs,run_data,not_strict):
    
 
    read_file_list = [read_file_name1,read_file_name2]
    
    split_files = []
    
    for file_ in read_file_list:
        
        split_files = split_files +split_file(file_,divide_by,run_dirs["dir_for_read_parts"],not_strict)
        
        
    #print split_files
    #print "==================\n"
  
    final_file_list = []
    done_final_files = []
    for split_f1 in split_files:
                    
        for split_f2 in split_files:
            
            if split_f1 == split_f2:
                continue
            
            tmp_ck = 0
            for final_files in done_final_files:
                if split_f1[-2:] == final_files[-2:] and split_f1[:3] == final_files[:3]:
                    
                    tmp_ck = 1
                   
            
                            
            if split_f1[-2:] == split_f2[-2:] and tmp_ck == 0:
                final_file_list.append([os.path.join(run_dirs["dir_for_read_parts"],split_f1),os.path.join(run_dirs["dir_for_read_parts"],split_f2)])
                done_final_files.append(split_f1)
                
    return final_file_list

def sort_output_files(file_list):
    # sort file chunks returned from the cluster according to the numbers after PART_
    
    dict_ = {}
    for file_ in file_list:
        
        dict_[int(file_.split("PART_")[1].split("_")[0])] = file_
        
    dic_list = dict_.keys()
    dic_list.sort()
    
    final_list = []

    for key in dic_list:
        
        print key
        final_list.append(dict_[key])
        
    return final_list

def cat_files(final_file_name,work_directory_master,file_count):
	
	#print final_file_name
	#print file_count
	#print work_directory_master
	#raw_input()
        final_file = open(final_file_name,"w")
	
	if len(file_count) > 1:
		file_count=sort_output_files(file_count)
    
        for bowtie_file in file_count:
            #print bowtie_file
            tmp_list = []
            
            tmp_count = 0
            for line in open(os.path.join(work_directory_master,bowtie_file)):
                
                if tmp_count == 1000000:
                    final_file.writelines(tmp_list)
                    tmp_list = []
                    tmp_count = 0
                  
                tmp_list.append(line)
                tmp_count += 1
            final_file.writelines(tmp_list)
            #print "Done",bowtie_file

def cat_files_sam(final_file_name,work_directory_master,file_count):
	
	#print final_file_name
	#print file_count
	#print work_directory_master
	#raw_input()
        final_file = open(final_file_name,"w")
	
	if len(file_count) > 1:
		file_count=sort_output_files(file_count)
	
	# only for @ headers    
        for bowtie_file in file_count:
            #print bowtie_file
            tmp_list = []
            
            tmp_count = 0
            for line in open(os.path.join(work_directory_master,bowtie_file)):
                if line.startswith("@"):
			
			if tmp_count == 1000000:
			    final_file.writelines(tmp_list)
			    tmp_list = []
			    tmp_count = 0
			  
			tmp_list.append(line)
			tmp_count += 1
		
		else:
			final_file.writelines(tmp_list)
			break
		
	for bowtie_file in file_count:
            #print bowtie_file
            tmp_list = []
            
            tmp_count = 0
            for line in open(os.path.join(work_directory_master,bowtie_file)):
                if line.startswith("@"):
			continue
		if tmp_count == 1000000:
		    final_file.writelines(tmp_list)
		    tmp_list = []
		    tmp_count = 0
		  
		tmp_list.append(line)
		tmp_count += 1
	    final_file.writelines(tmp_list)
	
	final_file.close()

def print_files(final_file_name,work_directory_master,file_count):
        cnt = 0
        meep = []
        
        file_count.sort()
    
        for bowtie_file in file_count:
            tmp_list = []
            
            for line in open(os.path.join(work_directory_master,bowtie_file)):
                
                print line
               
                
def get_queue_list(username="emmett"):

    cluster_proc = subprocess.Popen(['qstat','-u',username],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    q_list = []
    for line in cluster_proc.stdout:
        if line.startswith("job") or line.startswith("--"):
            continue
        
        q_list.append(line.strip().split())
    
    return q_list

def queue_info(total_files,tag,username="emmett"):
    
    q_list = get_queue_list(username)
    q_time_dict = {}
    
    q_jobs = []
    
    #take only jobs for this current script - identified by tag
    for q in q_list:
        if q[2].startswith(tag):
            
            q_jobs.append(q)
    
    total_jobs = len(q_jobs)
    
    names_jobs_present = []
    number_jobs = []
    for q1 in q_jobs:
        
        names_jobs_present.append(q1[2][:10])
        
        number_jobs.append(int(q1[2].split("_")[1]))
    
    numbers_not_present = []
    for i in range(11,total_files+11):
        if i not in number_jobs:
            numbers_not_present.append(i)
    
    return total_jobs,numbers_not_present,names_jobs_present
 
def time_now():
	curr_time = datetime.now()
	return curr_time.strftime("%c")
 
def create_dir(dir_,overwrite=1):
    
	if os.path.exists(dir_):
		if overwrite == 1:
			print >> sys.stderr,"WARNING: %s exists! Overwriting..."%(dir_)
			shutil.rmtree(dir_)
		else:
			return
	
	os.mkdir(dir_)

def build_cluster_Script(node_aligner_part,node_unmapped_part,aligner_part,unmapped_part,read_file_name1,read_file_name2,script_name,dir_job_done,files_created,local_copy,tmp_node_dir,diag=0):
	
			
	### Build bash script for each node!!!        
	output_files = []
	
	two_refs = False
	# create a temporary directory for use
	if diag != 1:
		script_cmd="if [ -d %s ]; then\n\tx=1\nelse\n\tmkdir %s\nfi\n"%(run_data["tmp_node_dir"],run_data["tmp_node_dir"])
	else:
		script_cmd=""
	# only engage local copy script if it is not in testmode
	if diag != 1:
		if local_copy == True:
			
			local_read1= os.path.join(tmp_node_dir,"local_read1")
			
			files_created.append(local_read1)
			
			if paired == True:
				local_read2= os.path.join(tmp_node_dir,"local_read2")
				
				script_cmd+="cp %s %s\n"%(read_file_name2,local_read2)
				files_created.append(local_read2)
				read_file_name2=local_read2
				
			script_cmd+="cp %s %s\n"%(read_file_name1,local_read1)
			
			read_file_name1=local_read1
		
	
			
	if run_data["read_pp"] != "none":
		read_pp = run_data["read_pp"]
		
		read_iter = [1]
		if paired == True:
			read_iter.append(2)
		
		
		for read_num in read_iter:
			
			readpp_step_number = len(run_data["read_pp"].strip().split("|"))
			print "steps: ",readpp_step_number
			
			
			pp_Read_prev_output = ""
				
				
			pp_read_prefix = os.path.join(os.path.split(node_aligner_part)[0],"read_processing_")
				
			read_step_number = 1
			
			# process reads into a format accepted by the aligner
			for line in run_data["read_pp"].strip().split("|"):
			    
				#first iteration use read name
				if read_step_number == 1:
					print "MOOP"
			    
				
					if read_num == 1:
						line =line.replace("pp_Read",read_file_name1)
					else:
						line =line.replace("pp_Read",read_file_name2)
			    
				      
				elif line.find("pp_In"):
					# remaining iterations use previous step output (pp_IN) as input 
					print line.find("pp_In")
					print line 
					if line.find("pp_In"):
						print "IN FOUND"
					   
						line =line.replace("pp_In",pp_Read_prev_output)
						print line
						print "end"
				else:
					print >> sys.stderr,"ERROR: No input (pp_In) marker found in this command:\n%s\nPlease check your config file."%(line)
					sys.exit()
					    
					
				    # for final step make a final file                
				if readpp_step_number == read_step_number:
					final_read_name = "%sfinal_Read%s"%(pp_read_prefix,str(read_num))
					
					if read_num == 1:
						read_file_name1 = final_read_name
					else:
						read_file_name2 = final_read_name
							    
					line =line.replace("pp_Out",final_read_name)
				       
					files_created.append(final_read_name)
				       
				       
				    # temporary output is created which is used as input    
				elif line.find("pp_Out"):
					tmp_read_name = "%sRead%s_Part%s"%(pp_read_prefix,str(read_num),str(read_step_number))
					line =line.replace("pp_Out",tmp_read_name)
					files_created.append(tmp_read_name)
				       
					
					pp_Read_prev_output = tmp_read_name
				else:
					print >> sys.stderr,"ERROR: No output (pp_Out) marker found in this command:\n%s\nPlease check your config file."%(line)
					sys.exit()
				    
				script_cmd+= line+"\n"     
				read_step_number+=1
	if diag == 1:
		script_cmd+="touch %s\n"%(os.path.join(dir_job_done,"1"))
	"""
	if run_data["ref_pp"] != "none":
		print >> sys.stderr,"No reference processing selected."
			
	else:
	"""
	final_ref1 = aligner_data[aligner]["ref_file1"]
	try:
		final_ref2 = aligner_data[aligner]["ref_file2"]
		two_refs = True
		#print >> sys.stderr,"Second reference found."
	except:
		#print >> sys.stderr,"WARNING: No second reference found."
		do = "nothing"
	
	if diag == 1:
		script_cmd+="touch %s\n"%(os.path.join(dir_job_done,"2"))
		
	analysis_str = ""
	if paired == True:
		print "Paired-end analysis"
		
		analysis_str = "paired_node_analysis"
	else:
		
		analysis_str = "single_node_analysis"

	# add analysis steps to node script
	# replace all markers with specific files, depending on the number of analysis steps
	analysis_steps = len(run_data[analysis_str].strip().split("|"))
	step_number = 1
	prev_output1 = ""
	prev_output2 = ""
	prev_output3 = ""

	for line in run_data[analysis_str].strip().split("|"):
	    
		line =line.replace("ALIGNER_Ref1",final_ref1)
		if two_refs == True:
			
			line =line.replace("ALIGNER_Ref2",final_ref2)
		line =line.replace("ALIGNER_in1",prev_output1)
		#print "ALIGNER_in1",prev_output1
		
		if line.find("ALIGNER_in2") > -1:
			#print "ALIGNER_in2",prev_output2
			line =line.replace("ALIGNER_in2",prev_output2)
			
		if line.find("ALIGNER_in3") > -1:
			line =line.replace("ALIGNER_in3",prev_output3)
		
		
		#replace read markers with actual reads
		line =line.replace("ALIGNER_Read1",read_file_name1)
		
		if paired == True:
			line =line.replace("ALIGNER_Read2",read_file_name2)
		
		
		# for final step
		if step_number == analysis_steps:
			if line.find("ALIGNER_output1") > -1:
				out1 = "%s_final_1"%(node_aligner_part)
				line =line.replace("ALIGNER_output1",out1)
				files_created.append(out1)
				output_files.append(out1)
			
			elif prev_output1 != "":
				output_files.append(prev_output1)
			
			
			if line.find("ALIGNER_output2") > -1:
				out2 = "%s_final_2"%(node_aligner_part)
				line =line.replace("ALIGNER_output2",out2)
				files_created.append(out2)
				output_files.append(out2)
				
			elif prev_output2 != "":
				output_files.append(prev_output2)
				
			if line.find("ALIGNER_output3") > -1:
				out3 = "%s_final_3"%(node_aligner_part)
				line =line.replace("ALIGNER_output3",out3)
				files_created.append(out3)
				output_files.append(out3)
				
			elif prev_output3 != "":
				output_files.append(prev_output3)
	    
					    
		else:
			if line.find("ALIGNER_output1") > -1:
				out1 = "%s_%s_1"%(node_aligner_part,str(step_number))
				line =line.replace("ALIGNER_output1",out1)
				files_created.append(out1)
				prev_output1 = out1
			
			if line.find("ALIGNER_output2") > -1:
				out2 = "%s_%s_2"%(node_aligner_part,str(step_number))
				line =line.replace("ALIGNER_output2",out2)
				files_created.append(out2)
				prev_output2 = out2
								
			if line.find("ALIGNER_output3") > -1:
				out3 = "%s_%s_3"%(node_aligner_part,str(step_number))
				line =line.replace("ALIGNER_output3",out3)
				files_created.append(out3)
				prev_output3 = out3
				
							
			#line =line.replace("ALIGNER_output1",node_aligner_part+str(step_number))
			#line =line.replace("ALIGNER_unmapped",node_unmapped_part+str(step_number))
			#files_created.append(node_aligner_part+str(step_number))
			#files_created.append(node_unmapped_part+str(step_number))
			
			
			
			
		    
		script_cmd+= line+"\n"        
		step_number +=1
	
	
	#copy local to lustre mount
	if diag != 1:
		for cp_file in output_files:
			#print cp_file
			#print type(cp_file)
			#print cp_file[1:]
			#p
			#print "cp %s %s\n"%(cp_file,"moo")
			#print aligner_part
			#print type(aligner_part)
			#rint "%s_%s"%(aligner_part,str(cp_file)[-1])
			script_cmd+="cp %s %s\n"%(cp_file,"%s_%s"%(aligner_part,cp_file[-1]))  
		
	#clean up!
	for file_ in files_created:
		script_cmd+="rm %s\n"%(file_)
	#script_cmd+="rm %s\n"%(node_bowtie_part)
	
	#raw_input("RETURNED::\n"+script_cmd) 
	#script_name,dir_job_done
	if diag == 1:
		script_cmd+="touch %s\n"%(os.path.join(dir_job_done,"3"))
	else:
		script_cmd+="touch %s\n"%(os.path.join(dir_job_done,script_name[:-3]))
	
	return script_cmd

# collect user information before making scripts
def get_info():
        clusteralign_root = ""
        while 1:        
                clusteralign_root=raw_input("Please give the directory where you have extracted clusteralign. If you are currently running ClusterAlign from the directory only type 'cwd':\n")
                
                if clusteralign_root == "cwd":
                    clusteralign_root = os.getcwd()
                    break
                else:
                    if os.path.isdir(clusteralign_root):
                        print >> sys.stderr,"Directory confirmed"
                        break
                    else:
                        print >> sys.stderr,"Directory does not exist, please retry"
        master_tmp = ""
        while 1:        
                master_tmp=raw_input("Please provide a temporary directory on the master server where ClusterAlign can write.\nPLEASE NOTE: Large volumes of data will be temporarily stored here, ideally this should be on a LUSTRE MOUNT.\n Alternatively please consult your cluster administrator to determine the best location.\n\n Default temp directory is 'tmp' folder in the clusteralign directory)\nLeave blank to accept default:\n")
                                
                if master_tmp.strip() == "":
                        
                        master_tmp = os.path.join(clusteralign_root,"tmp")
                        break
                
                if os.path.isdir(master_tmp):
                        out = 0
                        while 1:
                                ans=raw_input("Directory exists, you wish to overwrite all data ?")
                                
                                if ans.lower() == "y":
                                        out = 1
                                        break
                                                                        
                                elif ans.lower() == "n":
                                        break
                        if out ==1:
                                break
        node_tmp = ""                             
                
        node_tmp=raw_input("Please provide a temporary directory on the CLUSTER NODE where ClusterAlign can write.\nPLEASE NOTE:Consult your cluster administrator to determine the best location.\nDEFAULT IS SET TO '/tmp/clusteralign':\n")
                         
        if node_tmp.strip() == "":
                node_tmp ="/tmp/clusteralign"
                
        
        username =raw_input("Please enter your cluster username:\n")
      
        cluster_q=raw_input("Please enter the cluster queue you want to be set as default:\n")
        
        mem = 0
        
        while 1:
            node_mem =raw_input("Please enter the node memory you want to be set as default (default value = 3gb)\nPlease use either mb or gb for units\nExample: 150mb\n")
            
            if node_mem.lower().endswith("mb"):
                mem = int(float(node_mem[:-2])*1048576)
                break
                
            elif node_mem.lower().endswith("gb"):
                
                mem = int(float(node_mem[:-2])*1073741824)
                break
            
            elif node_mem.strip() == "":
                
                mem = 3273741824
                break
            
            else:
                
                print >> sys.stderr,"Incorrect value entered."
    
        return clusteralign_root,username,cluster_q,mem,master_tmp,node_tmp


# crate scripts using script directory
def create_scripts():
        
        
        clusteralign_root,username,cluster_q,mem,master_tmp,node_tmp=get_info()

        current_dir = os.listdir(clusteralign_root)
        
        
        for file_ in current_dir:
            
                if file_.endswith(".config"):
                    
                        ans=raw_input("Please note: ALL default config files will be restored!! Do you want this ? Y/N\n")
                        
                        if ans.lower() == "y":
                                break
                                
                        
                        elif ans.lower() == "n":
                                sys.exit("Creating Config Files aborted...")
                                return
                
        
        scripts_dir =  os.path.join(clusteralign_root,"backup_scripts")
        software_dir = os.path.join(clusteralign_root,"aligners")
        test_dir = os.path.join(clusteralign_root,"test_data")
        
        
        config_files = os.listdir(scripts_dir)
        
        for config_file in config_files:
                
                write_list = []   
                
                for line in open(os.path.join(scripts_dir,config_file)):
                        
                        line = line.replace("input_CLUSTERQUEUE",cluster_q)
                        line = line.replace("input_USERNAME",username)
                        line = line.replace("input_NODEMEM","%i"%(mem))
                        line = line.replace("SOFTWAREDIR",software_dir)
                        line = line.replace("TESTDIR",test_dir)
                        
                        line = line.replace("input_MASTER",master_tmp)  
                        line = line.replace("input_NODE",node_tmp)
			line = line.replace("input_ROOT",clusteralign_root)
                        
                        
                        write_list.append(line)
                        
                open(config_file,"w").writelines(write_list)

# Read command line arguements
#################################################
###
##PARAMTERZ
######################################################################################    
log_list =[]



curr_time = datetime.now()
start_time = curr_time.strftime("%c")
log_list.append("Start time: %s\n"%(start_time))

start = ztime.time()

bowtie_line= []

nextline=0
afterargs = 0
unmapped = ""
unmap_switch = 0


read1_switch = 0
read2_switch = 0

std_out = False
paired = False

diag = False
aligner = None
local_copy = False
paired = False
input_out = None
input_length = None
not_strict = False
input_read1 = None
input_read2 = None
input_queue = None
make_script = None


if len(sys.argv) == 1:
    
	print >> sys.stderr,"ClusterAlign v0.15\n==================="
	print >> sys.stderr,"First time running this? Run '--makescripts' option to create in-built aligner scripts."
	print >> sys.stderr,"The ONLY command required to run ClusterAlign is the name of the aligners' config file. All essential information is present in this file. ie. for the built-in bowtie aligner.\nTHIS MUST BE THE FIRST ARGUEMENT!!:"
	print >> sys.stderr,"python cluster_align.py bowtie.config\n"
	print >> sys.stderr,"The following commands MUST be specified on the command line and are not in the config file:"
	print >> sys.stderr,"\n--paired : Run a paired end analysis"
	print >> sys.stderr,"\n--localcopy : Copy reads to local nodes before alignment. Recommended if cluster doesnt have a Lustre mount"
	print >> sys.stderr,"\n--testmode : Run a basic diagnostic to determine whether your aligner file works"
	print >> sys.stderr,"\n--makescripts : Create all built-in aligner scripts"
	print >> sys.stderr,"\n--notstrict : Do not verify reads are fastq - this option not recommended unless you are using input other than FASTQ (which is not inherently supported)"
	
	print >> sys.stderr,"\n\nTo overwrite parameters in the config file on the command line use the following commands:"
	print >> sys.stderr,"\n--length=x : change the number of reads per file to x\n\n--out=/path/to/output_directory :specify a folder for your output\n\n--read1=/path/to/readfile.fastq : specify the first read file\n\n--read2=/path/to/readfile2.fastq specifiy the second pair read file\n\n--queue=xx change cluster queue to xx"
	
	sys.exit()

aligner = sys.argv[1]

if os.path.isfile(aligner):
	print >> sys.stderr,"Aligner file found."
else:
	if aligner.startswith("--make"):
		pass
	else:
		print >> sys.stderr,"Aligner file '%s' not found! Please verify path and file name!"%(aligner)
		sys.exit()
	
	
for arg_ in sys.argv[1:]:
	
	if arg_.endswith("config"):
		continue
 
		
	if  arg_ == "--localcopy":
				
		local_copy = True
	    
	elif arg_ == "--paired":
		paired = True
	
	elif arg_ == "--testmode":
		diag = True
	
	elif arg_.startswith("--makescripts"):
		create_scripts()
		sys.exit("Exiting...")
		
	elif arg_.lower() == "--notstrict":
		not_strict = True
		
	elif arg_.startswith("--length="):
		input_length=arg_.split("=")[1]
	
	elif arg_.startswith("--out="):
		input_out=arg_.split("=")[1]
		
	elif arg_.startswith("--read1="):
		input_read1=arg_.split("=")[1]
		
	elif arg_.startswith("--read2="):
		input_read2=arg_.split("=")[1]
	
	elif arg_.startswith("--queue="):
		input_queue=arg_.split("=")[1]
	else:
		print >> sys.stderr,"Unrecognised paramter: %s"%(arg_)
		sys.exit()


		
config_aligner = "%s"%(aligner)
config_general = os.path.join(os.path.split(aligner)[0],"cluster.config")

aligner = aligner.split(".")[0]

#test if files exist
if os.path.isfile(config_aligner):
	print >> sys.stderr,"File found: %s"%(config_aligner)
else:
	sys.exit("ERROR: Aligner configuration file has not been found: %s\nAborting...."%(config_aligner))
	
if os.path.isfile(config_general):
	
	print >> sys.stderr,"File found: %s"%(config_general)
else:
	sys.exit("ERROR: Cluster configuration file has not been found: %s\nAborting...."%(config_general))


config_paths = [config_general,config_aligner]
##############################################################################################################
# SCRIPTING START
##############################################################################################################        

def check_analysis_lines(analysis_str,analysis_type):
	
	vocab_dict = {"read":["pp_Out","pp_Read"],"node":["ALIGNER_output1","ALIGNER_Ref1","ALIGNER_Read1"],"post":["Results_1","ALIGNER_output1"]}
	
	for key in vocab_dict:
		if key in analysis_type:
			essential_params = vocab_dict[key]
	for op_ in essential_params:
		
		if op_ not in analysis_str:
			print >> sys.stderr,"WARNING: FOLLOWING OPERAND MISSING! : %s \nPlease revisit command line for: %s\n"%(op_,analysis_type)
			sys.exit()
	
	
	



def read_config_file(config_paths,aligner,paired,input_out,log_list,input_read1,input_read2):
	# Read in config file
	############################################################################
	
	aligner_data = {}
	aligner_data[aligner] = {}
	cluster_data = {}
	run_data = {}
	
	read_pp=""
	ref_pp=""
	single_node_analysis = ""
	paired_node_analysis=""
	post_analysis = ""
	#get config options from file
	cfg_line_number= 0
	for config_path in config_paths:
		for line in open(config_path):
			cfg_line_number+=1
			if line.startswith("#"):
				continue
			
			if line.strip() == "":
				continue
			
			if line.find("=") == -1:
				print >> sys.stderr,"WARNING: Line %s does not contain expected syntax.\n%s\n Please Place a '#' in front of the line or include an '=' to show value."%(cfg_line_number,line)
				sys.exit()
			    
			# ALIGNER DATA
			################# 
			
			
			
			if line.strip().startswith("read_directory"):
				#print "}}}{}{}{}{}{",line
				
				if input_read1 != None:
					aligner_data[aligner][line.strip().split("=")[0].strip()] = "none"
					continue
				
				if line.strip().split("=")[1].strip().lower() == "none":
					
					aligner_data[aligner][line.strip().split("=")[0].strip()] = "none"
					continue
				
				elif os.path.isdir(line.strip().split("=")[1].strip()):
					read_file_name1 = line.strip().split("=")[1].strip()
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
					continue
				else:
					
					print >> sys.stderr,"Read_directory does not exist or is not a direcory:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					sys.exit()	
			
			
			elif line.strip().startswith("test_genome_file"):
				
				if os.path.isfile(line.strip().split("=")[1].strip()):
					read_file_name1 = line.strip().split("=")[1].strip()
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
					continue
				else:
					
					print >> sys.stderr,"test_genome_file does not exist:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					sys.exit()
			
			elif line.strip().startswith("test_read_file1"):
				
				if os.path.isfile(line.strip().split("=")[1].strip()):
					read_file_name1 = line.strip().split("=")[1].strip()
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
					continue
				else:
					
					print >> sys.stderr,"test_read_file1 does not exist:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					sys.exit()
			elif line.strip().startswith("test_read_file2"):
				if os.path.isfile(line.strip().split("=")[1].strip()):
					read_file_name1 = line.strip().split("=")[1].strip()
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
					continue
				else:
					
					print >> sys.stderr,"test_read_file2 does not exist:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					sys.exit()
					
			elif line.strip().startswith("read_file1"):
				
				if input_read1 != None:
					
					aligner_data[aligner][line.strip().split("=")[0].strip()] = input_read1
					continue
			    
				if os.path.isfile(line.strip().split("=")[1].strip()):
					read_file_name1 = line.strip().split("=")[1].strip()
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
					continue
				else:
					
					print >> sys.stderr,"Read_file1 does not exist:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					sys.exit()
			
			    
			elif line.strip().startswith("read_file2"):
				
				if input_read2 != None:
					
					aligner_data[aligner][line.strip().split("=")[0].strip()] = input_read2
					paired = False
					continue
				
				
					
				if line.strip().split("=")[1].strip() == "none":
					
					if paired == True:
						print >> sys.stderr,"Read_file2 is not specified. Paired-end Analysis is disabled."
						paired = False
						continue
				
				elif input_read1 != None:
					print >> sys.stderr,"Read_file2 is not specified on the commandline! Paired-end Analysis is disabled."
					paired = False
					continue
					
				else:
					#paired = True
					if os.path.isfile(line.strip().split("=")[1].strip()):
					   
						read_file_name2 = line.strip().split("=")[1].strip()
						aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
						
						continue
					else:
						print >> sys.stderr,"read_file2 does not exist:%s\nPlease verify path."%(line.strip().split("=")[1].strip())
					
						sys.exit()
			
			elif line.strip().startswith("read_prefix"):
				aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("read_suffix"):
				aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
				
			elif line.strip().startswith("paired_prefix_length"):
				try:
					aligner_data[aligner][line.strip().split("=")[0].strip()] = int(line.strip().split("=")[1].strip())
				except:
					print >> sys.stderr,"Error on line %s reading paired_prefix_length:\t%s\nThis should be an integer greater than 0"%(cfg_line_number,line.strip().split("=")[1].strip())
					sys.exit()
					
			    
			elif line.strip().startswith("ref_file1"):
				aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
				
			elif line.strip().startswith("ref_file2"):
				aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("output_directory"):
				if line.strip().split("=")[1].strip() == "none":
					std_out = True
				if input_out == None:
					aligner_data[aligner][line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
				else:
					aligner_data[aligner][line.strip().split("=")[0].strip()] = input_out
				    
			
			# RUN DATA
			   
			elif line.strip().startswith("read_pp"):
				read_pp+= line.strip().split("=")[1].strip()+"|"
			#run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("number_of_rp_steps"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("ref_pp"):
				ref_pp+=line.strip().split("=")[1].strip()+"|"
			
			#run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			
			elif line.strip().startswith("single_node_analysis"):
				single_node_analysis+=  line.strip().split("=")[1].strip()+"|"
				#run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("post_analysis"):
				post_analysis+=  line.strip().split("=")[1].strip()+"|"
				#run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
				
				
			elif line.strip().startswith("paired_node_analysis"):
				paired_node_analysis+=  line.strip().split("=")[1].strip()+"|"
				#run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("master_tmp_dir"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
				
				#run_data["master_tmp_dir"]
				run_data["job_done"] = os.path.join(line.strip().split("=")[1].strip(),"job_done")
				run_data["dir_for_read_parts"] = os.path.join(line.strip().split("=")[1].strip(),"read_parts")
				run_data["tmp_results_dir"] = os.path.join(line.strip().split("=")[1].strip(),"tmp_results_dir")
				
			#elif line.strip().startswith("tmp_results_dir"):   
			#    run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("tmp_node_dir"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("length_parameter"):
				if input_length == None:
					try:   
						run_data[line.strip().split("=")[0].strip()] = int(line.strip().split("=")[1].strip())
					    
					except:
						print >> sys.stderr,"Error on line %s reading length_parameter:\t%s\nThis should be an integer greater than 0"%(cfg_line_number,line.strip().split("=")[1].strip())
						sys.exit()
				else:
					try:   
						#len_param = input_length.split("=")[1]
						run_data[line.strip().split("=")[0].strip()] = int(input_length)
					    
					except:
						print >> sys.stderr,"Error, length_parameter given is not an integer"
						sys.exit()
			
			elif line.strip().startswith("script_dir"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("cluster_queue"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("username"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			    
			elif line.strip().startswith("node_mem"):   
				run_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			# CLUSTER DATA
			#################         
			    
			elif line.strip().startswith("cluster_platform"):   
				cluster_data[line.strip().split("=")[1].strip()] = {}
			
			elif line.strip().startswith("cluster_submit"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("cluster_submit_out"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("cluster_submit_error"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("queue_param"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("mem_param"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("cluster_acct"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("cluster_qacct_us"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			elif line.strip().startswith("cluster_js"):   
				cluster_data[line.strip().split("=")[0].strip()] = line.strip().split("=")[1].strip()
			
			else:
				print >> sys.stderr,"Error line %s contains no recognised variables. Please review.\n %s"%(cfg_line_number,line)
				sys.exit()
		
	run_data["read_pp"] =read_pp[:-1]
	run_data["ref_pp"] =ref_pp[:-1]

	run_data["single_node_analysis"] = single_node_analysis[:-1]
	run_data["paired_node_analysis"] = paired_node_analysis[:-1]
	
	run_data["post_analysis"] = post_analysis[:-1]
	
	if read_pp[:-1] != "none":
		check_analysis_lines(read_pp[:-1],"read_pp")
	
	check_analysis_lines(single_node_analysis[:-1],"single_node_analysis")
	
	if paired_node_analysis[:-1].lower() != "none":
		check_analysis_lines(paired_node_analysis[:-1],"paired_node_analysis")	
	else:
		print >> sys.stderr,"Paired-end Analysis disabled: No command provided."
		paired = False
	
	if post_analysis[:-1] != "none":
		check_analysis_lines(post_analysis[:-1],"post_analysis")	
	
	print >> sys.stderr,"Basic Parameter check successful."
	pp= pprint.PrettyPrinter(indent=4)    
	
	
	'''
	print "aligner_data"
	pp.pprint(aligner_data)
	print "cluster_data"
	pp.pprint(cluster_data)
	print "run_data"
	pp.pprint(run_data)
	'''
	print >> sys.stderr,"Parameters selected:\nOutput Directory: %s\nNumber of reads per file: %s\n"%(aligner_data[aligner]["output_directory"],run_data["length_parameter"])
	log_list.append("Parameters selected:\nOutput Directory: %s\nNumber of reads per file: %s\n"%(aligner_data[aligner]["output_directory"],run_data["length_parameter"]))	 
	   
	#raw_input()
	print >> sys.stderr,"[%s] Config file successfully parsed."%(time_now())
	log_list.append("[%s] Config file successfully parsed.\n"%(time_now()))
	############################################################################
	return run_data,aligner_data,cluster_data,log_list

def create_all_directories(run_data,aligner_data):
	############ Create all directories and store all final directories in a dictionary
	run_dirs = {}
		
	# tag is a random sequence to define each run
	tag = ''.join([choice(string.letters) for i in range(3)])
	
	
	# create root directories if they dont exist DO NOT OVERWRITE IF THEY DO
	create_dir(run_data["master_tmp_dir"],overwrite=0)
	create_dir(run_data["dir_for_read_parts"],overwrite=0)
	create_dir(run_data["tmp_results_dir"],overwrite=0)
	create_dir(run_data["job_done"],overwrite=0)
	
	
	
	run_dirs["script_dir"] = os.path.join(aligner_data[aligner]["output_directory"],"scripts")
	run_dirs["dir_job_done"] = os.path.join(run_data["job_done"],tag)
	run_dirs["dir_for_read_parts"] = os.path.join(run_data["dir_for_read_parts"],tag)
	run_dirs["tmp_results_dir"] = os.path.join(run_data["tmp_results_dir"],tag)
	run_dirs["cluster_out_dir"] = os.path.join(aligner_data[aligner]["output_directory"],"cluster_out")
	run_dirs["cluster_err_dir"] = os.path.join(aligner_data[aligner]["output_directory"],"cluster_error")
	run_dirs["diag_dir"] = os.path.join(run_data["tmp_results_dir"],"diagnostic_test")
	
	'''
	script_dir = os.path.join(aligner_data[aligner]["output_directory"],"scripts")
	dir_job_done = os.path.join(run_data["job_done"],tag)
	dir_for_read_parts= os.path.join(run_data["dir_for_read_parts"],tag)
	tmp_results_dir= os.path.join(run_data["tmp_results_dir"],tag)
	
	cluster_out_dir = os.path.join(aligner_data[aligner]["output_directory"],"cluster_out")
	cluster_err_dir =os.path.join(aligner_data[aligner]["output_directory"],"cluster_error")
	'''
	# create temp read files
	# create output directories for raw cluster logs       
	# check the output directory and if it exists -> exit
	create_dir(aligner_data[aligner]["output_directory"],overwrite=0)
	create_dir(run_dirs["cluster_out_dir"])
	create_dir(run_dirs["cluster_err_dir"])
	
	
	# make dir for bash scripts
	create_dir(run_dirs["script_dir"])
	# make dir for split read files
	create_dir(run_dirs["dir_for_read_parts"])
	# tmp_results_dir
	create_dir(run_dirs["tmp_results_dir"])

	#diag_dir = os.path.join(run_data["tmp_results_dir"],"diagnostic_test")
	# create directory for diagnostics	
	create_dir(run_dirs["diag_dir"])
	# jobs done
	create_dir(run_dirs["dir_job_done"])
	    
	print >> sys.stderr,"[%s] Temporary directories successfully created."%(time_now())
	
	return run_dirs,tag

def launch_script(paired,aligner_data,run_data,run_dirs,tag,local_copy,not_strict,log_list):
	
	########paramters
	read_file_name1 = aligner_data[aligner]["read_file1"]

	if paired == True:   
		read_file_name2 = aligner_data[aligner]["read_file2"]
	else:
		read_file_name2 = ""
	
	#######Split files
	############
	
	cluster_lines = []
	
	if paired == True:
		file_list = paired_split_files(read_file_name1,read_file_name2,run_data["length_parameter"],run_dirs,run_data,not_strict)
	else:
		file_list = split_file(read_file_name1,run_data["length_parameter"],run_dirs["dir_for_read_parts"],not_strict)
	
	#print file_list
	#raw_input()
		
	print >> sys.stderr,"[%s] Read files split into fragments of size %s."%(time_now(),run_data["length_parameter"])
	log_list.append("[%s] Read files split into fragments of size %s.\n"%(time_now(),run_data["length_parameter"]))
	
	
	total_files = len(file_list)
	#print total_files
	#print file_list
	
	unmapped_files = []
	count_unmap = 11
	bash_scripts = []
	job_list = []
	job_name_list = []
	
	## create a script for all file parts and submit to cluster
	for file_ in file_list:
	    
		#files included
		files_created = []
		
		# create file paths for unmapped files  
		if paired == True:
			node_unmapped_part=os.path.join(run_data["tmp_node_dir"],unmapped+str(count_unmap))
			node_aligner_part = os.path.join(run_data["tmp_node_dir"],os.path.split(file_[0])[1]+aligner)
							
			aligner_part = os.path.join(run_dirs["tmp_results_dir"],"%s_%s_%s"%(os.path.split(file_[0])[1],tag,aligner))
			unmapped_part =os.path.join(run_dirs["tmp_results_dir"],unmapped+str(count_unmap))
			unmapped_files.append(unmapped+str(count_unmap))
			
			script_name = "%s_%s_%s.sh"%(tag,count_unmap,aligner)
			script_cmd=build_cluster_Script(node_aligner_part,node_unmapped_part,aligner_part,unmapped_part,file_[0],file_[1],script_name,run_dirs["dir_job_done"],files_created,local_copy,run_data["tmp_node_dir"])
		

		else:
			node_unmapped_part=os.path.join(run_data["tmp_node_dir"],unmapped+str(count_unmap))
			node_aligner_part = os.path.join(run_data["tmp_node_dir"],os.path.split(file_)[1]+aligner)
							
			aligner_part = os.path.join(run_dirs["tmp_results_dir"],"%s_%s_%s"%(os.path.split(file_)[1],tag,aligner))
			unmapped_part =os.path.join(run_dirs["tmp_results_dir"],unmapped+str(count_unmap))
			unmapped_files.append(unmapped+str(count_unmap))

			script_name = "%s_%s_%s.sh"%(tag,count_unmap,aligner)
			script_cmd=build_cluster_Script(node_aligner_part,node_unmapped_part,aligner_part,unmapped_part,file_,read_file_name2,script_name,run_dirs["dir_job_done"],files_created,local_copy,run_data["tmp_node_dir"])
			
		#print script_cmd
		#raw_input()
		
		#save bash script
		bowtie_bash_script = os.path.join(run_dirs["script_dir"],script_name)
		bash_scripts.append(bowtie_bash_script)
		
		
		open(bowtie_bash_script,"w").write(script_cmd)
		#log_list.append(script_cmd)
		
		
		#### create out and error files for each job
		run_dirs["cluster_out_dir"],run_dirs["cluster_err_dir"]
		cluster_error_file = os.path.join(run_dirs["cluster_err_dir"],"%s_%s_%s.err"%(tag,count_unmap,aligner)) 
		cluster_output_file =os.path.join(run_dirs["cluster_out_dir"],"%s_%s_%s.out"%(tag,count_unmap,aligner)) 
		
		#print cluster_data["cluster_submit"],"-o","/home/hep/emmett/scripts/cluster/current_script/out.txt","-e","/home/hep/emmett/scripts/cluster/current_script/error.txt",cluster_data["queue_param"],run_data["cluster_queue"],cluster_data["mem_param"],'mem_free=%s'%(run_data["node_mem"]),bowtie_bash_script
		
		cluster_proc = subprocess.Popen([cluster_data["cluster_submit"],"-o",cluster_output_file,"-e",cluster_error_file,cluster_data["queue_param"],run_data["cluster_queue"],cluster_data["mem_param"],'mem_free=%s'%(run_data["node_mem"]),bowtie_bash_script],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		cluster_lines.append([cluster_data["cluster_submit"],"-o",cluster_output_file,"-e",cluster_error_file,cluster_data["queue_param"],run_data["cluster_queue"],cluster_data["mem_param"],'mem_free=%s'%(run_data["node_mem"]),bowtie_bash_script])
		
		for line in cluster_proc.stdout:
			job_list.append(line.split()[2])
			
			job_name_list.append(line.split('"')[1])
			#print >> sys.stderr,line 
		#print "out"
		
		
		#for line in cluster_proc.stderr:
		#    print line
		#print "ERRout"
		
		count_unmap+=1
		
	print >> sys.stderr,"[%s] Jobs submitted successfully"%(time_now())
	log_list.append("[%s] Jobs submitted successfully.\n"%(time_now()))
	return total_files,job_name_list,job_list,bash_scripts,cluster_lines,log_list

def check_file_sam(work_directory_master,file_count):
	not_sam = False
	for file_ in file_count:
		
		for line in open(os.path.join(work_directory_master,file_)):
			
			if line.strip() == "":
				continue
			
			
			if line.startswith("@"):
				break
			else:
				not_sam = True
				
		if not_sam == True:
			#print "not sam"
			return False
	#print "sam"
	return True

def manage_jobs(cluster_data,run_data,run_dirs,total_files,tag,job_name_list,bash_scripts,cluster_lines,log_list):

	progbar = progressBar(0,total_files)
	count_done_files = 0
	print >> sys.stderr,progbar
	add_to_stats_list = []
	# check to see when all jobs have completed on the cluster
	
	final_output_files = {}

	repeat_jobs = []
	while 1:
	    
		ztime.sleep(10)
		#print "erp"
		file_list = os.listdir(run_dirs["tmp_results_dir"])
		output_file_dict = {"1":[],"2":[],"3":[]}
		unmap_count = []
		
		progress_list = os.listdir(run_dirs["dir_job_done"])
		
		final_prog_list = []
		
		for prog in progress_list:
		
			final_prog_list.append(prog[:10])
			
			
		if len(final_prog_list) > count_done_files:
			progbar.updateAmount(len(final_prog_list))
			print >> sys.stderr,progbar
			 
		total_jobs,numbers_not_present,names_jobs_present = queue_info(total_files,tag)
		
		#print names_jobs_present
		#print final_prog_list
		#print job_name_list
		completed_jobs = []
		failed_jobs = []
		
		for job_ in job_name_list:
		
			if job_[:10] not in names_jobs_present:
			
				if job_[:10] in final_prog_list:
				
					completed_jobs.append(job_)
				
				else:
					failed_jobs.append(job_)
		
		if len(failed_jobs) > 0:
			print "FAIL",failed_jobs
			for fjob in failed_jobs:
				
				if fjob not in repeat_jobs:
					print >> sys.stderr,"This job has not completed successfully, resubmitting: %s\n"%(fjob)
					log_list.append("This job has not completed successfully, resubmitting: %s\n"%(fjob))
					#failed_jobs.append(fjob)
									
				
					resub_job = ""
					for cluster_list in cluster_lines:
						
						if fjob in cluster_list[-1]:
							resub_job = cluster_list
							print fjob
							print cluster_list
							
							break
			
					cluster_proc = subprocess.Popen(resub_job,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
								
								
								
								    
					for line in cluster_proc.stdout:
						#job_list.append(line.split()[2])
						job_name_list.append(line.split('"')[1])
						print >> sys.stderr,line
						
					repeat_jobs.append(fjob)
					continue		
								
				elif fjob in repeat_jobs:
					print >> sys.stderr,"This job has not completed successfully and resubmission has been unsuccessful, please review logs files: %s\n"%(fjob)
					log_list.append("This job has not completed successfully and resubmission has been unsuccessful, please review logs files: %s\n"%(fjob))
		
					
		if len(numbers_not_present) == total_files:
		
			if len(final_prog_list) == total_files:
				
				print >> sys.stderr,"[%s] Cluster data retrieved successfully."%(time_now())
				log_list.append("[%s] Cluster data retrieved successfully.\n"%(time_now()))
				for key in output_file_dict:
					for file_ in file_list:
				       
						if file_.endswith(key):
						    
							output_file_dict[key].append(os.path.join(run_dirs["tmp_results_dir"],file_))
						
											
				     
					# cat files into final alignment output file
					#print >> sys.stderr,"[%s] Creating unmapped files."%(time_now())
					#cat_files(os.path.join(aligner_data[aligner]["output_directory"],"umapped_reads."+aligner),run_dirs["dir_for_read_parts"],unmap_count)
				       
					#cat files into final unmapped output file
					print >> sys.stderr,"[%s] Creating output alignments."%(time_now())
					log_list.append("[%s] Creating output alignments.\n"%(time_now()))
					final_output_files[key] = os.path.join(aligner_data[aligner]["output_directory"],"output_%s.%s"%(aligner,key))
					
					#print output_file_dict[key]
					
					if check_file_sam(run_dirs["dir_for_read_parts"],output_file_dict[key]) == True:
						cat_files_sam(os.path.join(aligner_data[aligner]["output_directory"],"output_%s.%s"%(aligner,key)),run_dirs["dir_for_read_parts"],output_file_dict[key])
					else:
						cat_files(os.path.join(aligner_data[aligner]["output_directory"],"output_%s.%s"%(aligner,key)),run_dirs["dir_for_read_parts"],output_file_dict[key])
				
					
				break
		else:
			
			
			if len(failed_jobs) == total_files:
				
				print >> sys.stderr,"[%s] All jobs have failed. This indicates an error in the config script. Please inspect the error logs in %s.\n"%(time_now(),run_dirs["cluster_err_dir"])
				log_list.append("[%s] All jobs have failed. This indicates an error in the config script. Please inspect the error logs in %s.\n"%(time_now(),run_dirs["cluster_err_dir"]))
				break
			
			

	return add_to_stats_list,final_output_files,log_list


# read in config file

def testmode(run_dirs,aligner_data):
		
		error_lines = []
		output_lines = []
		script_name = "test.sh"
		done_files_diag = os.path.join(run_dirs["diag_dir"],"done_files")
		
		error_file= os.path.join(aligner_data[aligner]["output_directory"],"testmode.stderr")
		out_file = os.path.join(aligner_data[aligner]["output_directory"],"testmode.stdout")
		
		create_dir(done_files_diag)
		#script_cmd=build_cluster_Script(node_aligner_part,node_unmapped_part,aligner_part,unmapped_part,read_file_name1,read_file_name2,script_name,run_dirs["dir_job_done"],files_created,local_copy,run_data["tmp_node_dir"])
		script_cmd=build_cluster_Script(os.path.join(run_dirs["diag_dir"],"node_aligner_part"),os.path.join(run_dirs["diag_dir"],"node_aligner_part"),os.path.join(run_dirs["diag_dir"],"aligner_master_part"),os.path.join(run_dirs["diag_dir"],"aligner_master_part"),aligner_data[aligner]["test_read_file1"],aligner_data[aligner]["test_read_file2"],script_name,done_files_diag,[],local_copy,os.path.join(run_dirs["diag_dir"],"node_dir"),diag=1)
		
		
		
		test_script = os.path.join(run_dirs["diag_dir"],"test.sh")
		open(test_script,"w").writelines(script_cmd)
		
		diag_proc =subprocess.Popen(["bash",test_script],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		for line in diag_proc.stdout:
			error_lines.append(line)
			#print line
		
		for line in diag_proc.stderr:
			error_lines.append(line)
			#print line
		
		#raw_input("end")
		
		done_list = os.listdir(done_files_diag)
		print done_list
		while 1:	
			if "1" in done_list:
				print >> sys.stderr,"Read Processing successful."
			else:
				print >> sys.stderr,"Read Processing failed. Please inspect the following output logs:\n%s\n%s\n"%(error_file,out_file)
				break
		
			if "2" in done_list:
				print >> sys.stderr,"Reference Processing successful."
			else:
				print >> sys.stderr,"Reference Processing failed. Please inspect the following output logs:\n%s\n%s\n"%(error_file,out_file)
				break
				
			if "3" in done_list:
				print >> sys.stderr,"Analysis Processing successful."
				break
			else:
				print >> sys.stderr,"Analysis Processing failed. Please inspect the following output logs:\n%s\n%s\n"%(error_file,out_file)
				break
		
		final_output_files = {"1":"%s_final_1"%(os.path.join(run_dirs["diag_dir"],"node_aligner_part"))
		,"2":"%s_final_2"%(os.path.join(run_dirs["diag_dir"],"node_aligner_part"))
		,"3":"%s_final_3"%(os.path.join(run_dirs["diag_dir"],"node_aligner_part"))
		}
		
		launch_post_analysis_commands(run_data,run_dirs,aligner_data,final_output_files,paired)	
			
		open(error_file,"w").writelines(error_lines)
		open(out_file,"w").writelines(output_lines)
		
		shutil.rmtree(run_dirs["diag_dir"])
		ans=raw_input("Test mode completed. Continue? yes/no\n")
		
		if ans.lower() not in ["yes","y"]:
			sys.exit("Exiting....")
			
def arrange_read_files(aligner_data,paired):
	
	file_list = []
	paired_dict = {}
	for file_ in os.listdir(aligner_data[aligner]["read_directory"]):
		
		if file_.startswith(aligner_data[aligner]["read_prefix"]) and file_.endswith(aligner_data[aligner]["read_suffix"]):
			file_list.append(file_)
			
			if paired == True:
				pair_prefix = file_[:aligner_data[aligner]["paired_prefix_length"]]
				
				if pair_prefix not in paired_dict:
					paired_dict[pair_prefix] = [file_]
				else:
					paired_dict[pair_prefix].append(file_)
					
	if paired == True:
		
		for lane in paired_dict:
			if len(paired_dict[lane]) != 2:
				
				sys.exit("WARNING: Directory contains an odd number read without a pair. Directory should contain only single-end or paired-end reads!!! Offending pair is::\n%s\n Exiting...."%("\n".join(paired_dict[lane])))
		return file_list,paired_dict
	else:
		return file_list,paired_dict

"""	
final_ref = aligner_data[aligner]["ref_file"]
read_file_name1 = aligner_data[aligner]["read_file1"]
read_file_name1 = aligner_data[aligner]["read_file2"]

post_analysis_name = os.path.join(aligner_data[aligner]["output_directory"],"post_processing")


"""	
	
	
def build_post_analysis_script(run_data,aligner_data,final_output_files,diag,paired):
		
	post_script = ""
	
	if diag != True:
		final_ref = aligner_data[aligner]["ref_file1"]
		read_file_name1 = aligner_data[aligner]["read_file1"]
		
		if paired == True:
			read_file_name2 = aligner_data[aligner]["read_file2"]
		
		post_analysis_name = os.path.join(aligner_data[aligner]["output_directory"],"post_processing")
	else:
		final_ref = aligner_data[aligner]["test_genome_file"]
		read_file_name1 = aligner_data[aligner]["test_read_file1"] #os.path.join(run_dirs["diag_dir"],"node_aligner_part")#aligner_data[aligner]["read_file1"]
		read_file_name2 = aligner_data[aligner]["test_read_file2"] #os.path.join(run_dirs["diag_dir"],"node_aligner_part")#aligner_data[aligner]["read_file2"]
		
		post_analysis_name = os.path.join(run_dirs["diag_dir"],"post_processing")
	
	analysis_steps = len(run_data["post_analysis"].strip().split("|"))
	
	
	#print run_data["post_analysis"].strip().split("|")
	
	step_number = 1
	prev_output1 =""
	prev_output2=""
	prev_output3=""
	
	for line in run_data["post_analysis"].strip().split("|"):
	    
		line =line.replace("ALIGNER_Ref",final_ref)
		line =line.replace("ALIGNER_in1",prev_output1)
		
		if line.find("ALIGNER_in2") > -1:
			line =line.replace("ALIGNER_in2",prev_output2)
			
		if line.find("ALIGNER_in3") > -1:
			line =line.replace("ALIGNER_in3",prev_output3)
		
		### Input files from aligner	
		if line.find("Results_1") > -1:	
			#prev_output1 = final_output_files["1"] 	
			line =line.replace("Results_1",final_output_files["1"])
			
		if line.find("Results_2") > -1:
			#prev_output2 = final_output_files["2"]	
			line =line.replace("Results_2",final_output_files["2"])
			
		if line.find("Results_3") > -1:
			#prev_output3 = final_output_files["3"]
			line =line.replace("Results_3",final_output_files["3"])
			
			
			
		#replace read markers with actual reads
		line =line.replace("ALIGNER_Read1",read_file_name1)
		
		if paired == True:
			line =line.replace("ALIGNER_Read2",read_file_name2)
		
			
		# for final step
		if step_number == analysis_steps:
			out1 = "%s_final_1"%(post_analysis_name)
			line =line.replace("ALIGNER_output1",out1)
			#files_created.append(out1)
			#output_files.append(out1)
			
			if line.find("ALIGNER_output2") > -1:
				out2 = "%s_final_2"%(post_analysis_name)
				line =line.replace("ALIGNER_output2",out2)
				#files_created.append(out2)
				#output_files.append(out2)
				
			if line.find("ALIGNER_output3") > -1:
				out3 = "%s_final_3"%(post_analysis_name)
				line =line.replace("ALIGNER_output3",out3)
				#files_created.append(out3)
				#output_files.append(out3)
	    
					    
		else:
			if line.find("ALIGNER_output1") > -1:
				out1 = "%s_%s_1"%(post_analysis_name,str(step_number))
				line =line.replace("ALIGNER_output1",out1)
				prev_output1 = out1
				#files_created.append(out1)
			
			if line.find("ALIGNER_output2") > -1:
				out2 = "%s_%s_2"%(post_analysis_name,str(step_number))
				line =line.replace("ALIGNER_output2",out2)
				#files_created.append(out2)
				prev_output2 = out2
								
			if line.find("ALIGNER_output3") > -1:
				out3 = "%s_%s_3"%(post_analysis_name,str(step_number))
				line =line.replace("ALIGNER_output3",out3)
				#files_created.append(out3)
				prev_output3 = out3
	
	
		post_script+= line+"\n"        
		step_number +=1
	
	return post_script
		

def launch_post_analysis_commands(run_data,run_dirs,aligner_data,final_output_files,paired):
	if run_data["post_analysis"].strip().lower() == "none":
		print >> sys.stderr,"[%s] Post analysis processing disabled."%(time_now())
		return
	diag = ""
	print >> sys.stderr,"[%s] Performing post analysis Processing."%(time_now())
	
	post_script_cmd = build_post_analysis_script(run_data,aligner_data,final_output_files,diag,paired)
	
	#print post_script_cmd
	
	if diag != True:
		
		post_script = os.path.join(run_dirs["script_dir"],"post_analysis.sh")
		
	else:
		post_script = os.path.join(run_dirs["diag_dir"],"post_analysis.sh")
		
	error_lines=[]	
	open(post_script,"w").writelines(post_script_cmd)
	post_analysis_proc =subprocess.Popen(["bash",post_script],stderr=subprocess.PIPE)
	
	#for line in post_analysis_proc.stdout:
		#error_lines.append(line)
		#print line
	if diag == True:
		print >> sys.stderr,"Following output was recieved from post-analysis script test:\n"
		
	for line in post_analysis_proc.stderr:
		error_lines.append(line)
		
		if diag == True:
			print >> sys.stderr,line
	
	if diag == True:
		open(os.path.join(aligner_data[aligner]["output_directory"],"testmode_post_analysis.log"),"w").writelines(error_lines)
	else:
		open(os.path.join(aligner_data[aligner]["output_directory"],"post_analysis.log"),"w").writelines(error_lines)
	#raw_input("end")
		
def pipeline_instance(aligner_data,run_data,local_copy,paired,diag,log_list):
	# create all directories for analysis
	run_dirs,tag =create_all_directories(run_data,aligner_data)
	
	if diag == 1:
		
		testmode(run_dirs,aligner_data)
		diag = 0
	

	# Split reads, create a script for each read, launch scripts and return detail on launched jobs
	total_files,job_name_list,job_list,bash_scripts,cluster_lines,log_list = launch_script(paired,aligner_data,run_data,run_dirs,tag,local_copy,not_strict,log_list)

	# manage and resubmit failed jobs
	add_to_stats_list,final_output_files,log_list = manage_jobs(cluster_data,run_data,run_dirs,total_files,tag,job_name_list,bash_scripts,cluster_lines,log_list)
	
	#########  Statistics and Times
	finish_time = datetime.now()
	end = ztime.time()
	duration = finish_time - start_time
	dur_secs = end-start
	print >> sys.stderr, "Pipeline exectution complete [%s elapsed]"%(duration)
	print >> sys.stderr, "Seconds: [%s elapsed]"%(dur_secs)
	
	log_list.append("Seconds: [%s elapsed]"%(dur_secs))
	log_list.append("Pipeline exectution complete [%s elapsed]\n"%(duration))
	open(os.path.join(aligner_data[aligner]["output_directory"],tag+"_logfile_"+aligner+"_cluster.txt"),"w").writelines(log_list)  
	
	# get statistics for jobs run on cluster
	collect_stats(job_list,cluster_data,run_data,tag)
	
	launch_post_analysis_commands(run_data,run_dirs,aligner_data,final_output_files,paired)
	
	# delete all temporary dirs
	delete_tmp_dirs(run_dirs,run_data,diag)
	
	return final_output_files


################## main script part
#parse config file
run_data,aligner_data,cluster_data,log_list = read_config_file(config_paths,aligner,paired,input_out,log_list,input_read1,input_read2)

if input_queue != None:
	print >> sys.stderr,"Cluster Queue Changed to: %s"%(input_queue)
	run_data["cluster_queue"] = input_queue

start_time = datetime.now()        

log_list.append("length_parameter %s\n"%(run_data["length_parameter"]))


master_output_dir = aligner_data[aligner]["output_directory"]

create_dir(master_output_dir)



if aligner_data[aligner]["read_directory"].lower().strip() != "none":
	
	file_list,paired_dict = arrange_read_files(aligner_data,paired)	
	
	print "filelist",file_list
	print "paired",paired_dict
	#raw_input()

	if paired == True:
		cnt =0
		for lane in paired_dict:
			
			aligner_data[aligner]["output_directory"] = os.path.join(master_output_dir,lane)
			aligner_data[aligner]["read_file1"] = os.path.join(aligner_data[aligner]["read_directory"],paired_dict[lane][0])
			aligner_data[aligner]["read_file2"] = os.path.join(aligner_data[aligner]["read_directory"],paired_dict[lane][1])
			
			print aligner_data[aligner]["output_directory"]
			print aligner_data[aligner]["read_file1"] 
			print aligner_data[aligner]["read_file2"]
			
			if cnt > 0:
				diag=0
			final_output_files = pipeline_instance(aligner_data,run_data,local_copy,paired,diag,log_list)
			cnt+=1
	else:
		cnt=0
		for read_file in file_list:
			
			aligner_data[aligner]["output_directory"] = os.path.join(master_output_dir,os.path.split(read_file)[-1].split(".")[0])
			aligner_data[aligner]["read_file1"] = os.path.join(aligner_data[aligner]["read_directory"],read_file)
			aligner_data[aligner]["read_file2"] = "none"
			
			if cnt > 0:
				diag=0
			final_output_files = pipeline_instance(aligner_data,run_data,local_copy,paired,diag,log_list)
			cnt+=1

else:
	final_output_files = pipeline_instance(aligner_data,run_data,local_copy,paired,diag,log_list)

	






