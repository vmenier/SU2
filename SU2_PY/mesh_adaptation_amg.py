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
	
	SU2.run.adaptation_amg (config, partitions)
	
#: def mesh_adaptation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
