#####packages needed
#sistema
from os import walk,path,mkdir,listdir,remove,rmdir,makedirs
import warnings

#data analyze
from pandas import read_csv, DataFrame,concat,to_numeric
from numpy import array, arange, log, sqrt,meshgrid, rot90,linspace,sort
from scipy import interpolate
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
import igor


def load_nanonis_sts(file_path,media):
    """
    Loads a Nanonis STM spectroscopy data file and extracts the numerical data into a Pandas DataFrame.

    Parameters:
    file_path (str): Path to the Nanonis .dat file.

    Returns:
    tuple: (metadata, dataframe) where:
        - metadata (dict) contains experimental parameters.
        - dataframe (pd.DataFrame) contains the spectroscopy data.
    """
    with open(file_path, "r", encoding="latin1") as file:
        lines = file.readlines()

    # Extract metadata before [DATA] section
    metadata = {}
    data_start_index = None

    for i, line in enumerate(lines):
        if line.strip() == "[DATA]":
            data_start_index = i + 1
            break
        else:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                key, value = parts
                metadata[key] = value
            elif len(parts) > 2:
                metadata[parts[0]] = parts[1:]

    # Extract column headers
    columns = lines[data_start_index].strip().split("\t")

    # Extract numerical data
    data_lines = lines[data_start_index + 1:]
    data_values = [line.strip().split("\t") for line in data_lines if line.strip()]

    # Convert to DataFrame
    df = DataFrame(data_values, columns=columns).apply(to_numeric, errors='coerce')

    df = df.sort_values(by=columns[0])
    # Resetando o índice para manter a numeração sequencial
    df = df.reset_index(drop=True)

    col_ida = [];col_volta = []
    col_media_volta = None
    for col in df.columns:
        if 'Current' in col:
            if '[bwd]' in col:
                if '[AVG]' in col:
                    col_media_volta = col
                else:
                    col_volta.append(col)
            else:
                if '[AVG]' in col:
                    col_media_ida = col
                else:
                    col_ida.append(col)           
    V_col_ida = [];V_col_volta = []
    V_col_media_volta = None
    for col in df.columns:
        if 'Bias' in col and 'cal' not in col:
            #print(col)
            if '[bwd]' in col:
                if '[AVG]' in col:
                    V_col_media_volta = col
                else:
                    V_col_volta.append(col)
            else:
                if '[AVG]' in col:
                    V_col_media_ida = col
                else:
                    V_col_ida.append(col)         
    if media == True:
        return [[df[V_col_media_ida],df[col_media_ida]],[df[V_col_media_volta],df[col_media_volta]]]
    else:
        ida = [];volta = []
        for i in range(len(col_ida)):
            ida.append([df[V_col_ida[i]],df[col_ida[i]]])
            if len(col_ida)==(len(col_volta)):
                volta.append([df[V_col_volta[i]],df[col_volta[i]]])
            else:
                volta.append([df[V_col_ida[i]],df[col_ida[i]]])
        return [ida,volta]

def get_name_file(path):
    ct =0
    for i in range(len(path)):
        letra = path[len(path)-1-i]
        if letra =='/' or letra == '\\':
            ct = len(path)-i
            return path[ct:]
    return path[ct:]

def load_multiple_files_nanonis(files_names,media,grid,nx=1,ny=1):
    Hyper_ida = {'V':[],'I':[],'name':[]};Hyper_volta = {'V':[],'I':[]}
    for names in files_names:
        if media:
            sts_ida,sts_volta = load_nanonis_sts(names,media)
            Hyper_ida['V'].append(sts_ida[0]);Hyper_ida['I'].append(sts_ida[1]*pow(10,9));
            Hyper_volta['V'].append(sts_volta[0]);Hyper_volta['I'].append(sts_volta[1]*pow(10,9))
            Hyper_ida['name'].append(get_name_file(names))
        else:
            sts_ida,sts_volta = load_nanonis_sts(names,media)
            
            for i in range(len(sts_ida)):
                Hyper_ida['name'].append(get_name_file(names))
                #name = names +str(i)
                Hyper_ida['V'].append(sts_ida[i][0]);Hyper_ida['I'].append(sts_ida[i][1]*pow(10,9))
                Hyper_volta['V'].append(sts_volta[i][0]);Hyper_volta['I'].append(sts_volta[i][1]*pow(10,9))      
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
