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


def load_multiple_grid_nanosfurf(files_names):
    data = []
    for names in files_names:
        data.append(nid_read(names,verbose=False))
    return data

def inverter(y,tipo ='trocar'):
    if tipo =='mult':
        return -y
    elif tipo=='trocar':
        ynew = []
        for i in range(len(y)):
            ynew.append(y[len(y)-1-i])
        return np.array(ynew)
def sts_nid(spec,media):
    I=spec.data.Spec.Forward["Tip Current"]
    V=spec.data.Spec.Forward["Tip voltage"]
    I2=spec.data.Spec.Backward["Tip Current"]
    V2=spec.data.Spec.Backward["Tip voltage"]
    if media:
        for i in range(len(I)):
            if i ==0:
                I_media=I[0];I2_media=I2[0]
                V_media=V[0];V2_media=V2[0]
            else:
                I_media=(I_media+I[i])/2;I2_media=(I2_media+I2[i])/2
                V_media=(V_media+V[i])/2;V2_media=(V2_media+V2[i])/2
                if I2_media[0]>=0:
                    I_v = inverter(I2_media)
                else:
                    I_v = I2_media
        return [[V_media,I_media*pow(10,9)],[V2_media,I_v*pow(10,9)]]
    else:
        V_ida = [];I_ida = [];V_volta = [];I_volta = []
        for i in range(len(I)):
            V_ida.append(V[i]);I_ida.append(I[i]*pow(10,9))
            if I2[i][0]>=0:
                I_v = inverter(I2[i])
            else:
                I_v = I2[i]
            V_volta.append(V2[i]);I_volta.append(I_v*pow(10,9))
        return [[V_ida,I_ida],[V_volta,I_volta]]


def get_name_file(path):
    ct =0
    for i in range(len(path)):
        letra = path[len(path)-1-i]
        if letra =='/' or letra == '\\':
            ct = len(path)-i
            return path[ct:]
    return path[ct:]


def get_grid_information(specs,names,media,grid,nx,ny):
    Hyper_ida = {'V':[],'I':[],'name':[]};Hyper_volta = {'V':[],'I':[]}
    stm_image_list = []
    for i in range(len(specs)):
        xmax,ymax,zmax =specs[i].param.Scan.range['Value']
        Z = specs[i].data.Image.Forward["Z-Axis"]*(pow(10,9))
        Z = Z-Z.min()
        xscan,yscan,stm_image = [np.linspace(0,xmax*pow(10,9),len(Z)),np.linspace(0,ymax*pow(10,9),len(Z[0])),Z]
        
        if media:
            sts_ida,sts_volta =sts_nid(specs[i],media)
            Hyper_ida['V'].append(sts_ida[0]); Hyper_ida['I'].append(sts_ida[1])
            Hyper_volta['V'].append(sts_volta[0]); Hyper_volta['I'].append(sts_volta[1])
            stm_image_list.append([xscan,yscan,stm_image])
            Hyper_ida['name'].append(get_name_file(names[i]))
        else:
            sts_ida,sts_volta =sts_nid(specs[i],media)
            
            for j in range(len(sts_ida[0])):
                Hyper_ida['V'].append(sts_ida[0][j]); Hyper_ida['I'].append(sts_ida[1][j])
                Hyper_volta['V'].append(sts_volta[0][j]); Hyper_volta['I'].append(sts_volta[1][j])
                Hyper_ida['name'].append(get_name_file(names[i]))
                stm_image_list.append([xscan,yscan,stm_image])

    
    if grid==None:
        return [stm_image_list,{"Hyper_ida":Hyper_ida,"Hyper_volta":Hyper_volta,'name':Hyper_ida['name']}]
    else:
        if nx!=1 and ny!=1:#media:
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
                if i%2==0:
                    #MV_ida.append(col_V_ida)
                    MI_ida.append(col_I_ida)
                    #MV_volta.append(col_V_volta); 
                    MI_volta.append(col_I_volta)
                else:
                    MI_ida.append(inverter(col_I_ida))
                    #MV_volta.append(col_V_volta); 
                    MI_volta.append(inverter(col_I_volta))                   
            #MHyper_ida = [Hyper_ida['V'],MI_ida] #pode ser que de problema no futuro se o eixo V nao for o mesmo para todos os dados
            #MHyper_volta = [Hyper_volta['V'],MI_volta] #pode ser que de problema no futuro se o eixo V nao for o mesmo para todos os dados
            return [stm_image_list,{"Hyper_ida":Hyper_ida,"Hyper_volta":Hyper_volta,'name':Hyper_ida['name'],'mapa_ida':MI_ida,'mapa_volta:':MI_volta}]
        else:
            MI_ida = []
            MI_volta =[]
            for i in range(len(Hyper_ida['I'])):
                MI_ida.append(Hyper_ida['I'][i])
                MI_volta.append(Hyper_volta['I'][i])
            return [stm_image_list,{"Hyper_ida":Hyper_ida,"Hyper_volta":Hyper_volta,'name':Hyper_ida['name'],'mapa_ida':MI_ida,'mapa_volta:':MI_volta}]


