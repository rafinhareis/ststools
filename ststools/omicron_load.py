#####packages needed
#sistema
from os import walk,path,mkdir,listdir,remove,rmdir,makedirs
import warnings

#data analyze
from pandas import read_csv, DataFrame,concat
from numpy import array, arange, log, sqrt,meshgrid, rot90,linspace,sort
from scipy import interpolate
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
import igor
# Get the path to the igor module
igor_path = path.dirname(igor.__file__)

#print(f"The igor package is located at: {igor_path}")
file_to_edit = path.join(igor_path, "binarywave.py")

# Read the content of the file and replace the deprecated usage
with open(file_to_edit, 'r', encoding='utf-8') as file:
    content = file.read()

# Replace 'np.complex' with 'np.complex128'
content = content.replace('_numpy.complex,', 'complex,')
# Write the modified content back to the file
with open(file_to_edit, 'w', encoding='utf-8') as file:
    file.write(content)
warnings.filterwarnings('ignore')
from igor.binarywave import load as load_ibw

#plot
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from .open_nid import *

def get_mapa_from_ibw(file):
    data_temp = load_ibw(file)
    data = data_temp['wave']
    ndim = data['wave_header']['nDim']
    sfA = data['wave_header']['sfA']
    sfB = data['wave_header']['sfB']
    matrix = data['wData']
    matrix = rot90(matrix)
    

    #dimensoes
    xtick = linspace(0,sfA[0]*ndim[0],ndim[0])*pow(10,9)
    ytick = linspace(0,sfA[1]*ndim[1],ndim[1])*pow(10,9)
    #bias voltage range
    V_list = []
    for line in data['note'].decode('utf-8').split('\r'):
        if ('Spectroscopy points per curve' in line)==True:
            for j in range(len(line.split())):
                if line.split()[j] =='=':
                    try:
                        V_list.append(float(line.split()[j+1][:-1]))
                    except ValueError:
                        pass
    V = linspace(V_list[1],V_list[2],ndim[2])
    if matrix[0][0][0]>0:
        return [xtick,ytick,V,-matrix*pow(10,9)]
    else:
        return [xtick,ytick,V,matrix*pow(10,9)]

###auxiliary funcions
def open_file(name):
    file = open(name, 'r', encoding='utf-8',
                errors='ignore')  # open the file.
    # creat a list that will be full filed with the file with a line as a list element.
    file_list = []
    for line in file:  # open a loop that cover the file.
        line = line.strip('\n')  # drop out all '\n' contained in every line.
        # change the spaces for a element of a list, ex: 'the energy is' --> ['the','energy','is'].
        line = line.split()
        file_list.append(line)  # add the line in the list file_list.
    file.close()  # close de file.
    return file_list

def ramp(file):
    for line in file:
        if len(line)>1:
            if line[0] == '#' and line[1] == 'Ramp':
                if line[2] == 'start':
                    v_start = round(float(line[len(line)-2] ),6)
                elif line[2] == 'end':
                    v_end = round(float(line[len(line)-2] ),6)
            elif ('Automatic' in line)== True and ('reversal' in line)==True:
                volta = line[len(line)-1]
    return [v_start,v_end,volta]

def txt_to_datafram(file):
    X = [];Y = []
    for line in file:
        if len(line)>=1:
            if line[0]!='#':
                try:
                    x,y = list(map(lambda x:float(x),line))   
                    X.append(x);Y.append(y)
                except ValueError:
                    pass
    df = DataFrame({'V':X,'I':Y})
    
    return df

def get_curve_number(path):
    if ('Spectroscopy' in path):
        for i in range(len(path)):
            n =len(path)-2-i
            if path[n]=='y' and path[n+1]=='-':
                ct = n+3
                break
        for i in range(ct,len(path)):
            if path[i]=='-':
                ct2 = i
                break
        return path[ct:ct2]
    else:
        return path

def get_name_file(path):
    ct =0
    for i in range(len(path)):
        letra = path[len(path)-1-i]
        if letra =='/' or letra == '\\':
            ct = len(path)-i
            return path[ct:]
    return path[ct:]


def load_multiple_files_omicron(files_names,grid,nx=1,ny=1):
    Hyper_ida = {'V':[],'I':[],'name':[]};Hyper_volta = {'V':[],'I':[]}
    MHyper_ida={'V':[],'I':[],'name':[]}
    ibw_file = False
    for names in files_names:
        if names.endswith('.asc'):
            pass
        elif names.endswith('.ibw'):
            ibw_file = True
            mapa = get_mapa_from_ibw(names)
            
            
            for i in range(len(mapa[3])):
                for j in range(len(mapa[3][0])):
                    MHyper_ida['V'].append(mapa[2])
                    MHyper_ida['I'].append(mapa[3][i][j])
                    MHyper_ida['name'].append(get_name_file(names))

            if mapa[2][0]>0:
                marker1 = True
                marker2 = False
                
            else:
                marker2 = True
                marker1 = False

        else:
            stm_files = open_file(names)
            v_start,v_end,volta =  ramp(stm_files)
            if volta == 'No':
                df = txt_to_datafram(stm_files)
                df = df.sort_values(by='V')
                # Resetando o índice para manter a numeração sequencial
                df = df.reset_index(drop=True)
                name_file = get_name_file(names)
                if v_start<v_end:
                    vmin = v_start;vmax = v_end
                else:
                    vmin = v_end;vmax = v_start
                Hyper_ida['V'].append(df['V']);Hyper_ida['I'].append(df['I']*pow(10,9));Hyper_ida['name'].append(name_file)
                Hyper_volta['V'].append(df['V']);Hyper_volta['I'].append(df['I']*pow(10,9))
            elif volta == 'Yes':
                df = txt_to_datafram(stm_files)
                
                for i in range(1,len(df)):
                    if df.iloc[i][0]>=v_end:
                        ct = i+1
                        break
                try:
                    df_ida = df.iloc[0:ct]
                    df_volta = df.iloc[ct:]
                except UnboundLocalError:
                    ct =int(len(df)/2)
                    df_ida = df.iloc[0:ct]
                    df_volta = df.iloc[ct:]
                df_ida = df_ida.sort_values('V')
                df_ida=df_ida.reset_index(drop=True)
                df_volta = df_volta.sort_values('V')
                df_volta=df_volta.reset_index(drop=True)
                #name_file = get_curve_number(names)
                name_file = get_name_file(names)
                if v_start<v_end:
                    vmin = v_start;vmax = v_end
                else:
                    vmin = v_end;vmax = v_start
                Hyper_ida['V'].append(df_ida['V']);Hyper_ida['I'].append(df_ida['I']*pow(10,9));Hyper_ida['name'].append(name_file)
                Hyper_volta['V'].append(df_volta['V']);Hyper_volta['I'].append(df_volta['I']*pow(10,9))
            else:
                print('erro to upload file: '+names)   
    if ibw_file:
        if marker1 == True and marker2==False:
            return [[None,mapa[0],mapa[1]],{"Hyper_ida":MHyper_ida,"Hyper_volta":MHyper_ida,'name':MHyper_ida['name'],'mapa_ida':mapa[3],'mapa_volta':-mapa[3]}]
        elif marker1 == False and marker2==True:
            return [[None,mapa[0],mapa[1]],{"Hyper_ida":MHyper_ida,"Hyper_volta":MHyper_ida,'name':MHyper_ida['name'],'mapa_ida':mapa[3],'mapa_volta':mapa[3]}]
        else:
            print('Problem to upload .ibw. files')
    else: 
        if grid==None:
            return [[None,None,None],{"Hyper_ida":Hyper_ida,"Hyper_volta":Hyper_volta,'name':Hyper_ida['name']}]
        else:
            MI_ida = [];MV_ida = []
            MI_volta =[];MV_volta = []
            for i in range(nx):
                col_I_ida = [];col_V_ida = []
                col_I_volta = [];col_V_volta = []
                for j in range(ny):
                    #print(ny*i+j)
                    try:
                        #col_V_ida.append(Hyper_ida['V'][ny*i+j])
                        col_I_ida.append(Hyper_ida['I'][ny*i+j])
                        #col_V_volta.append(Hyper_volta['V'][ny*i+j])
                        col_I_volta.append(Hyper_volta['I'][ny*i+j])                 
                    except IndexError:
                            print('Number of files are lower then elements of the grid rowsXcols.')
                MI_ida.append(col_I_ida)
                #MV_volta.append(col_V_volta); 
                MI_volta.append(col_I_volta)
                #MV_volta.append(col_V_volta); 
            return [[None,None,None],{"Hyper_ida":Hyper_ida,"Hyper_volta":Hyper_volta,'name':Hyper_ida['name'],'mapa_ida':MI_ida,'mapa_volta:':MI_volta}]