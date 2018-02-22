import os, sys
import _amgio as amgio
import numpy as np

def amg_call(config):
    
    cmd = '';
    cmd = "amg -in %s -sol %s -p 2 \
         -c %f -hgrad %.2f -hmin %le -hmax %le -out %s \
        -itp  %s  -nordg " \
        % (config['mesh_in'], config['sol_in'],  \
        config['size'],  config['hgrad'], config['hmin'], config['hmax'], \
        config['mesh_out'], config['itp_sol_in']);
        
    if config['adap_source'] != "":
        cmd += ' -source %s ' % config['adap_source']
    
    if config['adap_back'] != "":
        cmd += ' -back %s ' % config['adap_back']
    
    cmd += ' > %s' % config['amg_log']
    os.system(cmd);
    
def amg_call_python(mesh, config):
    
    remesh_options              = {}
    remesh_options['Lp']        = 2
    remesh_options['gradation'] = config['hgrad']
    remesh_options['target']    = config['size']
    remesh_options['logfile']   = config['amg_log']
    
    print mesh['sensor']
    
    ''' 
      TO ADD: 
     {'adap_back' 'hmax' 'hmin'
      'sol_in': 'current_sensor.solb', 'itp_sol_in': 'current.solb', 'metric_in': '', 'adap_source': '', 
     'mesh_in': 'current.meshb', 'mesh_out': 'current.new.meshb'}
    ''' 
    

# --- Read mesh using amgio module
def read_mesh(mesh_name, solution_name):
        
    Ver = [];
    Tri = [];
    Tet = [];
    Edg = [];
    Sol = [];
    SolTag = [];
    
    Markers = [];
    
    amgio.py_ReadMesh(mesh_name, solution_name, Ver, Tri, Tet, Edg, Sol, SolTag,  Markers);
        
    NbrTet = len(Tet)/5;
    Tet = np.reshape(Tet,(NbrTet, 5)).astype(int);

    NbrTri = len(Tri)/4;
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(int);
    
    NbrEdg = len(Edg)/3;
    Edg = np.reshape(Edg,(NbrEdg, 3)).astype(int);

    NbrVer = len(Ver)/3;
    Ver = np.reshape(Ver,(NbrVer, 3));
    
    SolSiz = len(Sol)/NbrVer;
    Sol = np.array(Sol).reshape(NbrVer,SolSiz).tolist();
    
    # First row of Markers contains dimension
    Dim = int(Markers[0]);
    
    mesh = dict();
    
    mesh['dimension']    = Dim
    
    if Dim == 3:
        mesh['xyz']          = Ver; 
    else :
        mesh['xy'] = np.stack((Ver[:,0],Ver[:,1]), axis=1)
    
        
    mesh['triangles']    = Tri;
    mesh['tetrahedra']   = Tet;
    mesh['edges']        = Edg;
    mesh['corners']      = [];
    mesh['solution']     = Sol;
    
    mesh['solution_tag'] = SolTag;
    
    mesh['id_solution_tag'] = dict();
    for i in range(len(SolTag)):
        mesh['id_solution_tag'][SolTag[i]] = i;
        
    mesh['markers']      = Markers;    
    
    return mesh;


def write_mesh(mesh_name, solution_name, mesh):
    
    Tri     = mesh['triangles']
    Tet     = mesh['tetrahedra']
    Edg     = mesh['edges']
    Sol     = mesh['solution']
    Markers = mesh['markers']
    
    Dim     = mesh['dimension']
    
    if "xyz" in mesh:
        Ver = mesh['xyz']
    else:
        NbrVer = len(mesh['xy']);
        z = np.transpose(np.zeros(NbrVer));
        Ver = np.c_[mesh['xy'],z]
    
    Ver = np.array(Ver).reshape(3*len(Ver)).tolist();
    Tri = np.array(Tri).reshape(4*len(Tri)).tolist();
    Tet = np.array(Tet).reshape(5*len(Tet)).tolist();
    Edg = np.array(Edg).reshape(3*len(Edg)).tolist();
    
    if len(Sol) > 1 :
        SolSiz = len(Sol[1])
        Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist();
    else:
        Sol = [];
    
    amgio.py_WriteMesh(mesh_name, solution_name, Ver, Tri, Tet, Edg, Sol, Markers, Dim);
    
    
def write_solution(solution_name, solution):
    
    Dim     = solution['dimension']
    Sol     = solution['solution']
    
    if "xyz" in solution:
        Ver = solution['xyz']
    else:
        NbrVer = len(solution['xy']);
        z = np.transpose(np.zeros(NbrVer));
        Ver = np.c_[solution['xy'],z]
    
    solution_tag = solution['solution_tag']
    
    NbrVer = len(Ver);
    Ver = np.array(Ver).reshape(3*len(Ver)).tolist();
    
    if len(Sol) > 1 :
        SolSiz = len(Sol[1])
        Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist();
    else:
        sys.stderr.write("## ERROR write_solution : No solution.\n");
        sys.exit(1);
        
    amgio.py_WriteSolution(solution_name, Ver, Sol, solution_tag, NbrVer, Dim)


def create_sensor(solution, sensor):
    
    if "xyz" in solution:
        Ver = solution['xyz']
    else:
        NbrVer = len(solution['xy']);
        z = np.transpose(np.zeros(NbrVer));
        Ver = np.c_[solution['xy'],z]
    
    #Ver = solution['xyz']; # Not needed for inria outputs
    
    NbrVer = len(Ver)
    Ver = np.array(Ver).reshape(3*len(Ver)).tolist();
    
    Dim = solution['dimension']
    Sol = np.array(solution['solution']);
    
    if sensor == "MACH":
        
        iMach = solution['id_solution_tag']['Mach'];
        sensor = Sol[:,iMach];
        sensor = np.array(sensor).reshape((len(sensor),1));
        sensor_header = ["Mach"];
        
    elif sensor == "PRES":
        
        iPres = solution['id_solution_tag']['Pressure'];
        sensor = Sol[:,iPres];
        sensor = np.array(sensor).reshape((len(sensor),1));        
        sensor_header = ["Pres"];
        
    elif sensor == "MACH_PRES":

        iPres  = solution['id_solution_tag']['Pressure'];
        iMach  = solution['id_solution_tag']['Mach'];
        mach   = np.array(Sol[:,iMach])
        pres   = np.array(Sol[:,iPres])
        sensor = np.stack((mach, pres), axis=1)    
        sensor_header = ["Mach", "Pres"];
                
    else :
        sys.stderr.write("## ERROR : Unknown sensor.\n");
        sys.exit(1);
    
    
    sensor_wrap = dict();
    
    sensor_wrap['solution_tag'] = sensor_header;
    
    if "xyz" in solution:
        sensor_wrap['xyz'] = solution['xyz'];
    else :
        sensor_wrap['xy'] = solution['xy'];
        
    sensor_wrap['dimension']    = solution['dimension'];
    sensor_wrap['solution']     = sensor;
    
    return sensor_wrap;


