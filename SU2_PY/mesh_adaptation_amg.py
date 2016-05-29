#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief Python script for running the anisotropic mesh adaptation using SU2 and AMG (Inria)
#  \author Victorien Menier
#  \version X.X.X
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import SU2
import os, time, sys, shutil, copy
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(): 

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=0,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--cycle", dest="cycle", default=1,
                      help="number of CYCLE adaptations", metavar="CYCLE")
    parser.add_option("-o", "--overwrite", dest="overwrite", default="False",
                      help="OVERWRITE_MESH the output mesh with the adapted one", metavar="OVERWRITE_MESH")
    parser.add_option("-s", "--save_all", dest="save_all", default="False",
                      help="SAVE_ALL the flow/adjoint/meshes solutions at each adaptation cycle", metavar="SAVE_ALL")

    (options, args)=parser.parse_args()

    options.partitions = int( options.partitions )
    options.cycle      = int( options.cycle      )
    options.overwrite  = options.overwrite == "True"    
    options.save_all   = options.save_all  == "True"
    
    # Run Mesh Adaptation
    mesh_adaptation ( options.filename   ,
                      options.partitions ,
                      options.cycle      ,
                      options.overwrite  ,
                      options.save_all    )

#: def main()


def amg_Inisol (config_adap, config):
	
	# TO DO HERE : DEAL WITH MESH FILE FORMAT
	#   If mesh format = inria and a solution is provided, continue
	#   If mesh format = SU2 and a solution is provided, call amg_Inisol to convert to the Inria format
	#   If a solution is not provided, call amg_Inisol which converts if needed
	
	rootDir    = config_adap.rootDir;
	adapDir    = config_adap.adapDir;
	partitions = config_adap.partitions;
	
	# --- Outputs: current mesh, restart, and sensor
	
	restart_name_dest  = "%s/current_restart.solb" % adapDir
	sensor_name_dest   = "%s/current_sensor.solb" % adapDir;  
	mesh_name_dest     = "%s/current.meshb" % adapDir;
	
	if config['ADAP_RESTART'] == 'NO':
		config_ini = copy.deepcopy(config);
				
		state = SU2.io.State()
		state.find_files(config_ini)

		# Set the number of partitions for parallel computations
		config_ini.NUMBER_PART = partitions
		config_ini.CONSOLE = 'QUIET'

		config_name = config_ini._filename;

		folder_inisol='initial_solution'; pull="../%s" % config_name; link="../%s" % config_ini['MESH_FILENAME'];

		config_ini.RESTART_FLOW_FILENAME = 'initial_sol.solb';
		config_ini.MESH_FORMAT = 'INRIA';
		config_ini.RESTART_SOL = 'NO';

		#config_ini.EXT_ITER = 1;
		
		jobNam = "SU2_ini.job"
		config_ini.OUTPUT_LOG = jobNam;
		print " Running SU2 on initial mesh \n Log: %s\n" % (jobNam) ;
				
		with SU2.io.redirect_folder(folder_inisol,pull,link):
			info = SU2.run.CFD(config_ini)
			shutil.copyfile (config_ini.RESTART_FLOW_FILENAME, restart_name_dest);  # restart solution
			shutil.copyfile ('mach.solb', sensor_name_dest );  # the solution used for error estimation (migh be different than the restart, ex: Mach)
			shutil.copyfile (config_ini.MESH_FILENAME        , mesh_name_dest   );
		print "...Done.\n\n";
	
	else :
		
		#--- An initial restart solution and an initial sensor are provided.
		#    -> Check their existence and copy them
		
		restart_name = "%s/%s" % (rootDir,config.ADAP_INI_RESTART_FILE);
		sensor_name  = "%s/%s" % (rootDir,config.ADAP_INI_SENSOR_FILE );
		mesh_name    = "%s/%s" % (rootDir,config.ADAP_INI_MESH_FILE   );
		
		error=0;
		
		if not os.path.isfile(restart_name):
			sys.stdout.write(' ## ERROR : No restart solution file was found.\n\n')
			error=1;
		
		if not os.path.isfile(sensor_name):
			sys.stdout.write(' ## ERROR : No restart sensor file was found.\n\n')
			error=1;
		
		if not os.path.isfile(mesh_name):
			sys.stdout.write(' ## ERROR : No restart mesh file was found.\n\n')
			error=1;
			
		if error == 1:
			sys.exit(1);
		
		shutil.copyfile (restart_name, restart_name_dest);   # restart solution
		shutil.copyfile (sensor_name , sensor_name_dest );   # the solution used for error estimation (migh be different than the restart, ex: Mach)
		shutil.copyfile (mesh_name   , mesh_name_dest   );
		
	
#: def amg_Inisol

def Parse_Adap_Options (config, config_adap):
	
	# --- Get the number of desired mesh adaptation loops
	#     and the mesh complexity for each of them
	
	config_adap.adap_complex = config['ADAP_COMPLEXITIES'].strip('()')
	config_adap.adap_complex = config_adap.adap_complex.split(",")
	
	config_adap.NbrIte = len(config_adap.adap_complex);
	
	if config_adap.NbrIte < 1 :
		print "  ## ERROR : Invalid number of iterations."
		sys.exit(1)
	
	# --- Get the number of sub-iterations for each mesh complexity
	
	config_adap.adap_subite = config['ADAP_SUBITE'].strip('()')
	config_adap.adap_subite = config_adap.adap_subite.split(",")
	
	if config_adap.NbrIte != len(config_adap.adap_subite) :
		print "  ## ERROR mesh adaptation: the number of mesh complexities (%d) and the number of sub-iterations (%d) don't match.\n" % (config_adap.NbrIte, config_adap.adap_subite)
		sys.exit(1)
	
	config_adap.NbrIteGlo = 0;
	for i in range(config_adap.NbrIte) :
		config_adap.adap_complex[i] = int(config_adap.adap_complex[i])
		config_adap.adap_subite[i]  = int(config_adap.adap_subite[i])
		config_adap.NbrIteGlo += config_adap.adap_subite[i];
	
	if 'ADAP_PATH' in config:
		path = "%s/" % config['ADAP_PATH'];
		os.environ['SU2_RUN'] = "%s:%s" % (os.environ['SU2_RUN'], path);
		os.environ['PATH'] = "%s:%s" % (os.environ['PATH'], path);
		config_adap.ADAP_PATH = path;
	else:
		config_adap.ADAP_PATH = "";
	
	if 'ADAP_BACK' in config:
		print " config['ADAP_BACK'] = %s\n" % config['ADAP_BACK'];
		if config['ADAP_BACK']=="YES":
			config_adap.ADAP_BACK = 1;
			
			if not 'ADAP_BACK_NAME' in config :
				print " ## ERROR : a back mesh name must be provided.\n ADAP_BACK option is ignored.\n";
				config_adap.ADAP_BACK = 0;
			else :
				config_adap.ADAP_BACK_NAME = config['ADAP_BACK_NAME'];
				
				#if not os.path.isfile(config_adap.ADAP_BACK_NAME):
				#	print " ## ERROR : Back mesh %s NOT FOUND!\n ADAP_BACK option is ignored.\n" %	config_adap.ADAP_BACK_NAME;
				#	config_adap.ADAP_BACK = 0;
		else:
			config_adap.ADAP_BACK = 0;
	else:
		config_adap.ADAP_BACK = 0;
			
	if (config_adap.ADAP_BACK):
		print "  -- Info : Use %s as a back mesh." % config_adap.ADAP_BACK_NAME;
	else:
		print "DONT GET BACK MESH";
	
	if 'ADAP_HIN' in config:
		config_adap.HMIN = float(config['ADAP_HMIN']);
	else:
		config_adap.HMIN = 0;
	
	if 'ADAP_HMAX' in config :
		config_adap.HMAX = float(config['ADAP_HMAX']);
	else:
		config_adap.HMAX = 1e300;
		
	if 'ADAP_HGRAD' in config:
		config_adap.HGRAD = float(config['ADAP_HGRAD']);
	else :
		config_adap.HGRAD = 3.0;
	
	
	# --- Print summary
	
	#print "Mesh adaptation summary: \n"
	#
	#for i in range(config_adap.NbrIte) :
	
# def: Parse_Adap_Options

def Call_AMG (config_adap, config):
	
	ite_glo = config_adap.ite_glo;
	ite_cpx = config_adap.ite_cpx;
	ite_sub = config_adap.ite_sub;
	
	outNam = "current.new.meshb";
	itpNam = "current_restart.itp.solb";
	
	rootDir    = config_adap.rootDir;
	
	cpx = config_adap.adap_complex[ite_cpx];
	
	hmin  = config_adap.HMIN;
	hmax  = config_adap.HMAX;
	hgrad = config_adap.HGRAD;
	
	jobNam = "AMG.%d.%.0f.job" % (ite_glo,cpx);
	
	if config_adap.ADAP_PATH != "":
		path = config_adap.ADAP_PATH;
	else:
		path = "";
	
	if os.path.isfile(outNam):
		os.remove(outNam);
		
	if os.path.isfile(itpNam):
		os.remove(itpNam);
		
	if (config_adap.ADAP_BACK):
		back = " -back %s/%s " % (rootDir, config_adap.ADAP_BACK_NAME);
	else:
		back = "";
	
	adapsrc = "%s/adap.source" % rootDir;
	if not os.path.exists(adapsrc):
		adapsrc = "./adap.source";
		open(adapsrc, 'a').close();
	
	amg_cmd = "%samg -in current.meshb -sol current_sensor.solb -p 2 -c %f -hgrad %.2f -hmin %le -hmax %le -out current.new.meshb -itp current_restart.solb %s -src %s > %s" % (path, cpx, hgrad, hmin, hmax, back, adapsrc, jobNam)
	
	print " Running AMG \n Log: %s\n" % (jobNam);
	
	os.system(amg_cmd);
	
	print "...Done.\n\n";
	
	if not os.path.isfile(outNam):
		msg = "  ## ERROR : AMG failed at global iteration %d (see %s)\n" % (ite_glo, jobNam);
		msg = "     Command : %s\n\n" % amg_cmd;
		sys.stdout.write(msg);
		sys.exit(1);
		
	if not os.path.isfile(itpNam):
		msg = "  ## ERROR AMG : solution interpolation failed  at global iteration %d (see %s)\n" % (ite_glo, jobNam);
		msg = "     Command : %s\n\n" % amg_cmd;
		sys.stdout.write(msg);
		sys.exit(1);
	
	os.rename(itpNam, "current.new_ini.solb");
	
	
	
#: def amg_RunIte

def Call_SU2 (config_adap, config):
	
	outNam = 'current.new_restart.solb';
	
	ite_glo = config_adap.ite_glo;
	ite_cpx = config_adap.ite_cpx;
	
	rootDir    = config_adap.rootDir;
	adapDir    = config_adap.adapDir;
	partitions = config_adap.partitions;
	
	
	cpx = config_adap.adap_complex[ite_cpx];
	
	config_cur = copy.deepcopy(config);

	state = SU2.io.State()
	state.find_files(config_cur)
	
	jobNam = "SU2.%d.%.0f.job" % (ite_glo,cpx);
	
	config_cur.NUMBER_PART = partitions
	config_cur.CONSOLE = 'QUIET'
	config_cur.OUTPUT_LOG = jobNam;
		
	config_cur.RESTART_SOL            = 'NO'
	config_cur.MESH_FILENAME          = 'current.new.meshb'
	config_cur.SOLUTION_FLOW_FILENAME = 'current.new_ini.solb';
	config_cur.RESTART_FLOW_FILENAME  = outNam;
	config_cur.MESH_FORMAT            = 'INRIA';

	#config_cur.EXT_ITER = 1;
	
	print " Running SU2 \n Log: %s\n" % (jobNam) ;
	
	if os.path.isfile(outNam):
		os.remove(outNam);
	
	info = SU2.run.CFD(config_cur)
	
	print "...Done.\n\n";
	
	if not os.path.isfile(outNam):
		msg = "  ## ERROR : SU2 failed at global iteration %d\n" % (ite_glo);
		sys.stdout.write(msg);
		sys.exit(1);
		
	shutil.copyfile ("mach.solb", "current.new_sensor.solb");
		
#: Call_SU2

# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation( filename             ,
                     partitions   = 0     , 
                     cycles       = 1     ,
                     overwrite    = False ,
                     save_all     = False  ):
								
	if not filename:
		print "  ## ERROR : a .cfg file must be provided."
		sys.exit(1)
	
	# --- Warning before deleting old files?	
	warn = 1;
		
	# ---------------------------------------------
	# --- Parse mesh adaptation options
	# ---------------------------------------------
	
	# ---  Set the name of the configuration file
	config_name = filename
	
	if not os.path.isfile(config_name):
		msg = "  ## ERROR : could not find configuration file %s\n\n" % config_name;
		sys.stdout.write(msg);
		sys.exit(1);
	
	# ---  Read the specified configuration file
	config = SU2.io.Config(config_name)
	
	config_adap = SU2.io.Config();
	Parse_Adap_Options(config, config_adap);
	
	config_adap.partitions = partitions;
	
	NbrIte       = config_adap.NbrIte;
	adap_complex = config_adap.adap_complex;
	adap_subite  = config_adap.adap_subite;
	
	# ---------------------------------------------
	# --- Create adaptation folder
	# ---------------------------------------------
	
	#for i in range(NbrIte) :
	#	print "Iteration %d : Complexity = %d, %d sub-iterations.\n" % (i,adap_complex[i], adap_subite[i])
		
	rootDir = os.getcwd();
	adapDir = "%s/ADAP" % rootDir;
         
	config_adap.rootDir = rootDir;
	config_adap.adapDir = adapDir;

	if os.path.exists(adapDir):
		sys.stdout.write('./ADAP exists. Removing old mesh adaptation in 10s.')
		sys.stdout.flush();
		if warn : time.sleep(10);
		shutil.rmtree(adapDir);
		sys.stdout.write(' Done!\n\n')
  
	os.makedirs(adapDir);
	os.chdir(adapDir);
	
	# ---------------------------------------------
	# --- Initial solution
	# ---------------------------------------------
	
	amg_Inisol (config_adap, config);
	# Out : current.meshb
	#       current_restart.solb
	#       current_sensor.solb
	
	# ---------------------------------------------
	# --- Run mesh adaptation loop
	# ---------------------------------------------
	
	ite_glo = 1;
	for ite_cpx in range(NbrIte):
				
		NbrSub = adap_subite [ite_cpx];
		Cpx    = adap_complex[ite_cpx];
		
		config_adap.ite_cpx = ite_cpx;
		
		plu=""
		if NbrSub>1 : plu="s"
		
		print " -- Mesh complexity %.0f: %d sub-iteration%s.\n" % (adap_complex[ite_cpx], NbrSub, plu)
		
		for ite_sub in range(NbrSub) :
			
			config_adap.ite_sub = ite_sub; 
			
			print "    -- Sub-iteration %d\n" % ite_sub
			
			config_adap.ite_glo = ite_glo;
			
			#--- Mesh adaptation and solution interpolation
			
			# In : current.meshb
			#      current_restart.solb
			#      current_sensor.solb
			Call_AMG(config_adap, config);
			# Out : current.new.meshb
			#       current.new_ini.solb
			
			#--- CFD computation
			
			# In : current.new.meshb
			#      current.new_ini.solb
			Call_SU2(config_adap, config);
			# Out : current.new_restart.solb
			#       current.new_sensor.solb
			
			
			#--- Save files
			shutil.copyfile ("current.new.meshb", "ite.%d.%.0f.meshb"%(ite_glo, Cpx));
			shutil.copyfile ("current.new_sensor.solb", "ite.%d.%.0f.solb"%(ite_glo, Cpx));
			shutil.copyfile ("current.new_restart.solb", "ite.%d.%.0f_restart.solb"%(ite_glo, Cpx));
			
			
			historyNam = "history.dat"
			if os.path.exists(historyNam):
				shutil.copyfile (historyNam, "ite.%d.%.0f_history.dat"%(ite_glo, Cpx));
				
			#--- Rename files for the next iteration
			shutil.move("current.new.meshb"       , "current.meshb");
			shutil.move("current.new_restart.solb", "current_restart.solb");
			shutil.move("current.new_sensor.solb" , "current_sensor.solb");		
			
			ite_glo = ite_glo+1;
						
		#: End sub-iteration
	
	#: ite_cpx
		
	sys.exit(1)
	
	
#: def mesh_adaptation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
