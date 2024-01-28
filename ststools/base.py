#####packages needed
from os import walk,path,mkdir,listdir,remove,rmdir
from pandas import read_csv, DataFrame,concat
from numpy import array, arange, log, sqrt,meshgrid, rot90,linspace,sort
from scipy import interpolate
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import ipywidgets as widgets
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter.messagebox import showinfo
from PIL import Image, ImageTk
from urllib.request import urlopen
import warnings
warnings.filterwarnings('ignore')
#global var
global folder
global filenames
global files_data
folder = None  # Declare global variable first
filenames = None
global hist_path
hist_path = 'not_path'
global hist_path_metal
hist_path_metal = 'not_path'
global hist_path_molecula
hist_path_molecula = 'not_path'
global curve_glob
global smooth_glob
global delta_glob
global selec 
curve_glob = 0; smooth_glob = 0; delta_glob = 0
global save_var
global stm
save_var = True
selec = {}

global save_path
global hist_path_forwad
global hist_path_metal_forwad
global hist_path_molecula_forwad
global hist_path_backward
global hist_path_metal_backward
global hist_path_molecula_backward

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
                    v_start = float(line[len(line)-2] )
                elif line[2] == 'end':
                    v_end = float(line[len(line)-2] )
            elif ('Automatic' in line)== True and ('reversal' in line)==True:
                volta = line[len(line)-1]
    return [v_start,v_end,volta]

def txt_to_datafram(file):
    X = [];Y = []
    for line in file:
        if len(line)>=1:
            if line[0]!='#':
                x,y = list(map(lambda x:float(x),line))   
                X.append(x);Y.append(y)
    return DataFrame({'V':X,'I':Y})

def get_curve_number(path):
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


def load_file(path,which='Nanosurf'):
    if which == 'Nanosurf':
        try:
            df = read_csv ( path, sep= ';', names= ['x','y','z'])
            first_value = df[df.columns[0]].iloc[0] 
            V = []
            for i in range(1,len(df)):
                if df[df.columns[0]].iloc[i]  == first_value:
                    count = i
                    #V.append(df[df.columns[0]].iloc[i])
                    break
                else:
                    V.append(df[df.columns[0]].iloc[i])

            n = int(len(df)/count)

            dataframes = []
            names = []
            Vmin = min(V)
            Vmax = max(V)
            for j in range(n):
                sts=[]
                volta =[]

                for i in range(len(df)):
                    if i<count :
                        filtro_alto = 18*pow(10,-9)
                        filtro_baixo = -filtro_alto
                        voltagem = df[df.columns[0]].iloc[i]
                        corrente = df[df.columns[2]].iloc[i+j*count]
                        if corrente>filtro_baixo and corrente < filtro_alto:
                            volta.append(voltagem)
                            sts.append(corrente)
                        else:
                            pass
                    else:
                        break

                V_new = volta
                try:
                    
                    if min(V_new)>Vmin:
                        Vmin = min(V_new)
                    if max(V_new)<Vmax:
                        Vmax = max(V_new)
                    df_new = DataFrame(V_new,columns = ['V'])
                    df_new['I'+str(j)] = sts
                    names.append(str(j) )
                    dataframes.append(df_new)
                except ValueError:
                    print('File %s with bad data.'%path)

            return {'ixv':[dataframes,names,[Vmin,Vmax]],'ixv_idaevolta':[[],[],[]]}
        except:
            df = read_csv ( path, sep= ';', names= ['x','y'])
            for i in range(len(path)):
                if path[i]=='/':
                    ct = i
                    break
            first_value = df[df.columns[0]].iloc[0] 
            V = array(df['x'])
            dataframes = []
            names = []
            Vmin = min(V)
            Vmax = max(V)
            for j in range(1):
                sts=[]
                volta =[]
                for i in range(len(df)):
                    if i<count :
                        filtro_alto = 18*pow(10,-9)
                        filtro_baixo = -filtro_alto
                        voltagem = df[df.columns[0]].iloc[i]
                        corrente = df[df.columns[2]].iloc[i+j*count]
                        if corrente>filtro_baixo and corrente < filtro_alto:
                            volta.append(voltagem)
                            sts.append(corrente)
                        else:
                            pass
                    else:
                        break

                V_new = volta
                try:
                    
                    if min(V_new)>Vmin:
                        Vmin = min(V_new)
                    if max(V_new)<Vmax:
                        Vmax = max(V_new)
                    df_new = DataFrame(V_new,columns = ['V'])
                    df_new['I'] = sts
                    names.append( path[ct+1:-4])
                    dataframes.append(df_new)
                except ValueError:
                    print('File %s with bad data.'%path)

            return {'ixv':[dataframes,names,[Vmin,Vmax]],'ixv_idaevolta':[[],[],[]]}
    elif which == 'Omicron':
        grid = path.endswith('.asc')
        if grid:
            pass

        else:
            stm_files = open_file(path)
            v_start,v_end,volta =  ramp(stm_files)
            if volta == 'No':
                df = txt_to_datafram(stm_files)
                name = get_curve_number(path)
                if v_start<v_end:
                    vmin = v_start;vmax = v_end
                else:
                    vmin = v_end;vmax = v_start
                return {'ixv':[[df],[name],[vmin,vmax]],'ixv_idaevolta':[[],[],[]]}
            elif volta == 'Yes':
                df = txt_to_datafram(stm_files)
                for i in range(1,len(df)):
                    if df.iloc[i][0]==v_end:
                        ct = i+1
                        break
                df_ida = df.iloc[0:ct]
                df_volta = df.iloc[ct:]
                df_volta = df_volta.sort_values('V')
                df_volta=df_volta.reset_index(drop=True)
                name = get_curve_number(path)
                if v_start<v_end:
                    vmin = v_start;vmax = v_end
                else:
                    vmin = v_end;vmax = v_start
                return {'ixv':[[df_ida],[name],[vmin,vmax]],'ixv_idaevolta':[[df_volta],[name],[vmin,vmax]]}    
            else:
                print('erro to upload file: '+path)   

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 80, fill = '█',init =0, final =100):
    """
    Call in a loop to create terminal progress bar
    Courtesy of S.O. user Greenstick, https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    
    percent = ("{0:." + str(decimals) + "f}").format(final * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
    	print()

def regression(x,y):
    x = array(x).reshape((-1, 1))
    model = LinearRegression().fit(x, y)
    b = model.intercept_
    a = model.coef_
    return [a[0],b]

def cut(x,y,xmin,xmax):
    x = array(x);y=array(y)
    xnew = [];ynew = []
    if x[0]>0:
        df = DataFrame({'x':x,'y':y})
        df = df.sort_values(by = 'x')
        x = array(df['x']);y = array(df['y'])
    for i in range(len(x)):
        if x[i]>=xmin and x[i]<=xmax:
            xnew.append(x[i])
            ynew.append(y[i])
    return [array(xnew),array(ynew)]

def integral(x,y):
    h = x[1]-x[0]
    I = (y[0]+y[len(y)-1])/2
    for i in range(1,len(y)-1):
        I+=y[i]
    I=h*I
    return I

def didv(x,y):
    h = x[1]-x[0]
    deri_y = []; deri_x =[]
    for i in range(1,len(x)-1):
        d = (y[i+1]-y[i-1])/(2*h)
        deri_x.append(x[i])
        deri_y.append(d)
    deri_y=array(deri_y)
    return [array(deri_x),deri_y/deri_y.max()]

def i_V(x,y):
  I =[]; V= []
  for i in range(1,len(x)-1):
    V.append(x[i])
    I.append(y[i])
  return [array(V),array(I)]

def gap_type(dx,dy,delta):
    f = interpolate.interp1d(dx, dy)
    marker1 = False; marker2 = False

    for i in range(len(dx)-1):
        if dx[0]<0:
            if dx[i]<=0 and dx[i+1]>=0:
               indice_x0 = i
               break
    marker1 = False; marker2 = False

    for i in range(len(dx)):
        if dx[0]<0:
          if i <= indice_x0:
            j= indice_x0-i
            if f(dx[j])<= delta and marker1 == False:
                xmin = dx[j]
            else:
               marker1 = True

          elif i > indice_x0:
            j= i 
            if f(dx[j])<= delta and marker2 == False:
                xmax = dx[j]
            else:
               marker2 = True
               
    try:
        gap = xmax-xmin
    except UnboundLocalError:
        xmin = 0
        xmax = 0
        gap = 0
    typ = abs(xmax) - abs(xmin)
    return [round(gap,2),round(typ,2),xmin,xmax]

def inter_x(x,y,z,dx = 100, dy= 100):

    interp = interpolate.RegularGridInterpolator((x,y),z)
    xnew = linspace(x.min(),x.max(),dx)
    ynew = linspace(y.min(),y.max(),dy)
    
    M_int = []
    for i in range(len(xnew)):
        pts = []
        for j in range(len(ynew)):
            pts.append([xnew[i],ynew[j]])
        pts = array(pts)
        col = interp(pts)
        M_int.append(array(col))
    
    return [xnew,ynew,array(M_int)]

def to_table(path_files,files_df,name = ''):
    #files = listdir(path_files)
    arq_I = open(path_files+name+'Dataframe_I_complete.csv','w')
    arq_didv = open(path_files+name+'Dataframe_dIdv_complete.csv','w')
    lines_I = []
    lines_didv = []
    marker=False
    marker2 = True
    for k in range(len(files_df)):
    #for item in files:
    #    if '.txt' in item:
            #df = read_csv(path_files+item)
            #number = ''
            #for car in item:
            #    if car=='_':
            #        break
            #    number+=car
            #col = df.columns
            df = files_df[k]
            number = str(k)
            col = df.columns
            x=df[col[0]];y=df[col[1]];dy=df[col[2]]
            dy=dy/dy.max()
            if marker == False:
                tam = len(x)
                xlimmin = x.min()
                xlimmax = x.max()
                lines_I.append([col[0]+'_'+number,',',col[1]+'_'+number])
                lines_didv.append([col[0]+'_'+number,',',col[2]+'_'+number])
                marker = True
                for j in range(len(x)):
                    lines_I.append([x[j],',',y[j]])
                    lines_didv.append([x[j],',',dy[j]])
            else:
                if x.min()>=xlimmin:
                    xlimmin=x.min()
                if x.max()<=xlimmax:
                    xlimmax=x.max()
                lines_I[0] = lines_I[0]+ [',',col[0]+'_'+number,',',col[1]+'_'+number]
                lines_didv[0] = lines_didv[0]+ [',',col[0]+'_'+number,',',col[2]+'_'+number]

                if len(x)<=tam:
                    for k in range(len(x)):
                        lines_I[k+1] = lines_I[k+1]+ [',',x[k],',',y[k]]
                        lines_didv[k+1] = lines_didv[k+1]+ [',',x[k],',',dy[k]]
                else:
                    for k in range(len(x)):
                        if k <tam:
                            lines_I[k+1] = lines_I[k+1]+ [',',x[k],',',y[k]]
                            lines_didv[k+1] = lines_didv[k+1]+ [',',x[k],',',dy[k]]
                        else:
                            lines_I.append([x[j],',',y[j]])
                            lines_didv.append([x[j],',',y[j]])
                    tam=len(x)
    for line in lines_I:
        line = list(map(lambda x:str(x),line))
        l = ''
        for c in line:
            l+=c 
        arq_I.writelines(l+'\n')
    arq_I.close()
    for line in lines_didv:
        line = list(map(lambda x:str(x),line))
        l = ''
        for c in line:
            l+=c 
        arq_didv.writelines(l+'\n')
    arq_didv.close()

    for k in range(len(files_df)):
    #for item in files_df:
        #if '.txt' in item:
            #df = read_csv(path_files+item)
            
            #number = ''
            #for car in item:
            #    if car=='_':
            #        break
            #    number+=car
            df = files_df[k]
            number = str(k)
            col = df.columns
            x=df[col[0]];y=df[col[1]];dy=df[col[2]]
            dy=dy/dy.max()
            if marker2:
                xnew  = arange(xlimmin,xlimmax,0.01)
                f = interpolate.interp1d(x,y)
                g = interpolate.interp1d(x,dy)
                df2 = DataFrame({'V':xnew,'I_'+number:f(xnew)})
                df3 = DataFrame({'V':xnew,'didv_'+number:g(xnew)})
                marker2=False
            else:
                f = interpolate.interp1d(x,y)
                g = interpolate.interp1d(x,dy)
                df2['I_'+number]=f(xnew)
                df3['didv_'+number]=g(xnew)
            printProgressBar(k,len(files_df[0])-1)
    df2 = df2.set_index('V')
    df3 = df3.set_index('V')
    df2.to_csv(path_files+name+'Dataframe_I_limeted_by_V.csv')
    df3.to_csv(path_files+name+'Dataframe_didv_limeted_by_V.csv')


def select_sts(path_file, file,n = [],smooth = 5,delta = 10,resolution = 0.01,save_um_file = False,metal = False,molecula = False,xmin_hist=-1,xmax_hist=1,name_table=''):
            global hist_path
            global hist_path_metal
            global hist_path_molecula
            global save_path
            #file = load_file(path_file)
            
            folder_save = 'sts_saves'
            for i in range(1,len(path_file[:-4])):
                if path_file[len(path_file[:-4]) -i] == '/':
                    ct = len(path_file[:-4]) -i
                    break
            try:
                #folder_name = path.join(folder_save,path_file[:-4][ct:])
                folder_name= folder_save+path_file[:-4][ct:]
            except UnboundLocalError:
                ct = 0
                #folder_name = path.join(folder_save,path_file[ct:-4])
                folder_name= folder_save+'/'+ path_file[:-4][ct:]
            try: 
                mkdir(folder_save )
            except FileExistsError:
                pass
                
            paste =folder_name
            save_path = paste
            try:
                mkdir(paste)
            except FileExistsError:
                for root, dirs, files_list in walk(paste):
                    for f in files_list:
                        #if file.endswith('.txt'):
                        remove(path.join(root, f))
                          
            hist ={'Curve':[],'Gap(V)':[] , 'Dop(Type)':[],'Dop(Value)':[]}
            files_df = []
            for i in range(len(file[0])):
                if (i in n) == False:
                        name = path.join(paste,file[1][i])
                        x= file[0][i][file[0][i].columns[0]];y = file[0][i][file[0][i].columns[1]]*pow(10,9)
                        p = int(smooth*len(y)/100)
                        if p%2==0:
                            p+=1
                            y = savgol_filter(y,p,1)
                        elif p==0:
                            pass
                        else:
                            y = savgol_filter(y,p,1)
                        dx,dy = didv(x,y)
                        gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                        hist['Curve'].append(i)
                        hist['Gap(V)'].append(round(gap,3))
                        hist['Dop(Value)'].append(round(typ,3))
                        if abs(round(typ,3))<=resolution:
                            tipo = 'neutro'
                        elif typ<-resolution:
                            tipo = 'n'
                        else:
                            tipo = 'p'
                        hist['Dop(Type)'].append(tipo)
                        f = interpolate.interp1d( file[0][i][file[0][i].columns[0]],file[0][i][file[0][i].columns[1]])
                        g = interpolate.interp1d(dx,dy)
                        xnew = arange(dx.min(),dx.max()-0.01,0.01)
                        df_new = DataFrame({'V':xnew,'I(nA)':f(xnew),'didv':g(xnew),'gap': str(gap),'tipo':str(tipo)}   )
                        files_df.append(df_new)
                        df_new = df_new.set_index('V')
                        if save_um_file:
                            if ct ==0:
                                df_new.to_csv(name +'_'+path_file[:-4][ct:]+'.txt') 
                            else:
                                df_new.to_csv(name +'_'+path_file[:-4][ct+1:]+'.txt')
            df_hist = DataFrame(hist)
            df_hist = df_hist.set_index('Curve')
            df_hist.to_csv(paste+'/'+'histogram.csv')
            hist_path = paste+'/'+'histogram.csv'
            if metal:
                    
                hist_metal = {'a':[],'a_neg':[],'a_pos':[]}
                dfs = file[0]
                for i in range(len(file[0])):                        
                        columns = dfs[i].columns
                        x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                        x_cut,y_cut = cut(x,y,xmin_hist,xmax_hist )
                        a,b = regression(x_cut,y_cut)
                        a = round(a,3); b = round(b,3)
                        hist_metal['a'].append(a)
                        dx,dy = didv(x,y)
                        x_cut_neg,y_cut_neg = cut(dx,dy,xmin_hist,0-resolution )
                        a_neg,b_neg = regression(x_cut_neg,y_cut_neg)
                        a_neg = round(a_neg,3); b_neg = round(b_neg,3)
                        x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist )
                        a_pos,b_pos = regression(x_cut_pos,y_cut_pos)
                        a_pos = round(a_pos,3); b_pos = round(b_pos,3)    
                        hist_metal['a_neg'].append(a_neg)
                        hist_metal['a_pos'].append(a_pos)
                df = DataFrame(hist_metal)
                df = df.to_csv(paste+'/'+'histogram_alphas.csv')
                hist_path_metal = paste+'/'+'histogram_alphas.csv'
            if molecula:
                his_mole = {'I_pos':[],'I_neg':[]}
                dfs = file[0]
                for i in range(len(file[0])):   
                    columns = dfs[i].columns
                    x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                    x_cut,y_cut = cut(x,y,xmin_hist ,xmax_hist  )
                    dx,dy = didv(x,y)
                    #I = round(integral(x_cut,y_cut),4)
                    x_cut_neg,y_cut_neg = cut(dx,dy,xmin_hist ,0-resolution )
                    I_neg = round(integral(x_cut_neg,y_cut_neg),4)
                    x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist  )           
                    I_pos = round(integral(x_cut_pos,y_cut_pos),4)
                    his_mole['I_pos'].append(I_pos)
                    his_mole['I_neg'].append(I_neg)
                df = DataFrame(his_mole)
                df = df.to_csv(paste+'/'+'histogram_molecule.csv')
                hist_path_molecula = paste+'/'+'histogram_molecule.csv'
            to_table(paste+'/',files_df,name=name_table)
            print("arquivos salvos na pasta "+ paste)

def select_sts_ida(path_file, file,n = [],smooth = 5,delta = 10,resolution = 0.01,save_um_file = False,metal = False,molecula = False,xmin_hist=-1,xmax_hist=1,name_table=''):
            global hist_path
            global hist_path_forwad
            global hist_path_metal_forwad
            global hist_path_molecula_forwad
            global save_path
            #file = load_file(path_file)

            folder_save = 'sts_saves'
            for i in range(1,len(path_file[:-4])):
                if path_file[len(path_file[:-4]) -i] == '/':
                    ct = len(path_file[:-4]) -i
                    break
            try:
                #folder_name = path.join(folder_save,path_file[:-4][ct:])
                folder_name= folder_save+path_file[:-4][ct:]
            except UnboundLocalError:
                ct = 0
                #folder_name = path.join(folder_save,path_file[ct:-4])
                folder_name= folder_save+'/'+ path_file[:-4][ct:]
            try: 
                mkdir(folder_save )
            except FileExistsError:
                pass
                
            paste =folder_name
            save_path = paste
            try:
                mkdir(paste)
            except FileExistsError:
                for root, dirs, files_list in walk(paste):
                    for f in files_list:
                        #if file.endswith('.txt'):
                        remove(path.join(root, f))
                          
            hist ={'Curve':[],'Gap(V)':[] , 'Dop(Type)':[],'Dop(Value)':[]}
            files_df = []
            for i in range(len(file[0])):
                if (i in n) == False:
                        name = path.join(paste,file[1][i])
                        x= file[0][i][file[0][i].columns[0]];y = file[0][i][file[0][i].columns[1]]*pow(10,9)
                        p = int(smooth*len(y)/100)
                        if p%2==0:
                            p+=1
                            y = savgol_filter(y,p,1)
                        elif p==0:
                            pass
                        else:
                            y = savgol_filter(y,p,1)
                        dx,dy = didv(x,y)
                        gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                        hist['Curve'].append(i)
                        hist['Gap(V)'].append(round(gap,3))
                        hist['Dop(Value)'].append(round(typ,3))
                        if abs(round(typ,3))<=resolution:
                            tipo = 'neutro'
                        elif typ<-resolution:
                            tipo = 'n'
                        else:
                            tipo = 'p'
                        hist['Dop(Type)'].append(tipo)
                        f = interpolate.interp1d( file[0][i][file[0][i].columns[0]],file[0][i][file[0][i].columns[1]])
                        g = interpolate.interp1d(dx,dy)
                        xnew = arange(dx.min(),dx.max()-0.01,0.01)
                        df_new = DataFrame({'V':xnew,'I(nA)':f(xnew),'didv':g(xnew),'gap': str(gap),'tipo':str(tipo)}   )
                        files_df.append(df_new)
                        df_new = df_new.set_index('V')
                        if save_um_file:
                            if ct ==0:
                                df_new.to_csv(name +'_'+'F_'+path_file[:-4][ct:]+'.txt') 
                            else:
                                df_new.to_csv(name +'_'+'F_'+path_file[:-4][ct+1:]+'.txt')
            df_hist = DataFrame(hist)
            df_hist = df_hist.set_index('Curve')
            df_hist.to_csv(paste+'/'+'histogram_forward.csv')
            hist_path_forwad = paste+'/'+'histogram_forward.csv'
            if metal:
                    
                hist_metal = {'a':[],'a_neg':[],'a_pos':[]}
                dfs = file[0]
                for i in range(len(file[0])):                        
                        columns = dfs[i].columns
                        x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                        x_cut,y_cut = cut(x,y,xmin_hist,xmax_hist )
                        a,b = regression(x_cut,y_cut)
                        a = round(a,3); b = round(b,3)
                        hist_metal['a'].append(a)
                        dx,dy = didv(x,y)
                        x_cut_neg,y_cut_neg = cut(dx,dy,xmin_hist,0-resolution )
                        a_neg,b_neg = regression(x_cut_neg,y_cut_neg)
                        a_neg = round(a_neg,3); b_neg = round(b_neg,3)
                        x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist )
                        a_pos,b_pos = regression(x_cut_pos,y_cut_pos)
                        a_pos = round(a_pos,3); b_pos = round(b_pos,3)    
                        hist_metal['a_neg'].append(a_neg)
                        hist_metal['a_pos'].append(a_pos)
                df = DataFrame(hist_metal)
                df = df.to_csv(paste+'/'+'histogram_alphas_forward.csv')
                hist_path_metal_forwad = paste+'/'+'histogram_alphas_forward.csv'
            if molecula:
                his_mole = {'I_pos':[],'I_neg':[]}
                dfs = file[0]
                for i in range(len(file[0])):   
                    columns = dfs[i].columns
                    x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                    dx,dy = didv(x,y)
                    x_cut,y_cut = cut(x,y,xmin_hist ,xmax_hist  )
                    #I = round(integral(x_cut,y_cut),4)
                    x_cut_neg,y_cut_neg = cut(dx,dy,xmin_hist ,0-resolution )
                    I_neg = round(integral(x_cut_neg,y_cut_neg),4)
                    x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist  )           
                    I_pos = round(integral(x_cut_pos,y_cut_pos),4)
                    his_mole['I_pos'].append(I_pos)
                    his_mole['I_neg'].append(I_neg)
                df = DataFrame(his_mole)
                df = df.to_csv(paste+'/'+'histogram_molecule_forward.csv')
                hist_path_molecula_forwad = paste+'/'+'histogram_molecule_forward.csv'
            to_table(paste+'/',files_df,name=name_table)
            print("arquivos salvos na pasta "+ paste)


def select_sts_volta(path_file, file,n = [],smooth = 5,delta = 10,resolution = 0.01,save_um_file = False,metal = False,molecula = False,xmin_hist=-1,xmax_hist=1,name_table = ''):
            global hist_path_backward
            global hist_path_metal_backward
            global hist_path_molecula_backward
            #file = load_file(path_file)
            
            folder_save = 'sts_saves'
            for i in range(1,len(path_file[:-4])):
                if path_file[len(path_file[:-4]) -i] == '/':
                    ct = len(path_file[:-4]) -i
                    break
            try:
                #folder_name = path.join(folder_save,path_file[:-4][ct:])
                folder_name= folder_save+path_file[:-4][ct:]
            except UnboundLocalError:
                ct = 0
                #folder_name = path.join(folder_save,path_file[ct:-4])
                folder_name= folder_save+'/'+ path_file[:-4][ct:]
            try: 
                mkdir(folder_save )
            except FileExistsError:
                pass
                
            paste =folder_name
            try:
                mkdir(paste)
            except FileExistsError:
                pass
                          
            hist ={'Curve':[],'Gap(V)':[] , 'Dop(Type)':[],'Dop(Value)':[]}
            files_df = []
            for i in range(len(file[0])):
                if (i in n) == False:
                        name = path.join(paste,file[1][i])
                        x= file[0][i][file[0][i].columns[0]];y = file[0][i][file[0][i].columns[1]]*pow(10,9)
                        p = int(smooth*len(y)/100)
                        if p%2==0:
                            p+=1
                            y = savgol_filter(y,p,1)
                        elif p==0:
                            pass
                        else:
                            y = savgol_filter(y,p,1)
                        dx,dy = didv(x,y)
                        gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                        hist['Curve'].append(i)
                        hist['Gap(V)'].append(round(gap,3))
                        hist['Dop(Value)'].append(round(typ,3))
                        if abs(round(typ,3))<=resolution:
                            tipo = 'neutro'
                        elif typ<-resolution:
                            tipo = 'n'
                        else:
                            tipo = 'p'
                        hist['Dop(Type)'].append(tipo)
                        f = interpolate.interp1d( file[0][i][file[0][i].columns[0]],file[0][i][file[0][i].columns[1]])
                        g = interpolate.interp1d(dx,dy)
                        xnew = arange(dx.min(),dx.max()-0.01,0.01)
                        df_new = DataFrame({'V':xnew,'I(nA)':f(xnew),'didv':g(xnew),'gap': str(gap),'tipo':str(tipo)}   )
                        files_df.append(df_new)
                        df_new = df_new.set_index('V')
                        if save_um_file:
                            if ct ==0:
                                df_new.to_csv(name +'_'+'B_'+path_file[:-4][ct:]+'.txt') 
                            else:
                                df_new.to_csv(name +'_'+'B_'+path_file[:-4][ct+1:]+'.txt')
            df_hist = DataFrame(hist)
            df_hist = df_hist.set_index('Curve')
            df_hist.to_csv(paste+'/'+'histogram_backward.csv')
            hist_path_backward = paste+'/'+'histogram_backward.csv'
            if metal:
                    
                hist_metal = {'a':[],'a_neg':[],'a_pos':[]}
                dfs = file[0]
                for i in range(len(file[0])):                        
                        columns = dfs[i].columns
                        x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                        x_cut,y_cut = cut(x,y,xmin_hist,xmax_hist )
                        a,b = regression(x_cut,y_cut)
                        a = round(a,3); b = round(b,3)
                        hist_metal['a'].append(a)
                        dx,dy = didv(x,y)
                        x_cut_neg,y_cut_neg = cut(dx,dy,xmin_hist,0-resolution )
                        a_neg,b_neg = regression(x_cut_neg,y_cut_neg)
                        a_neg = round(a_neg,3); b_neg = round(b_neg,3)
                        x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist )
                        a_pos,b_pos = regression(x_cut_pos,y_cut_pos)
                        a_pos = round(a_pos,3); b_pos = round(b_pos,3)    
                        hist_metal['a_neg'].append(a_neg)
                        hist_metal['a_pos'].append(a_pos)
                df = DataFrame(hist_metal)
                df = df.to_csv(paste+'/'+'histogram_alphas_backward.csv')
                hist_path_metal_backward = paste+'/'+'histogram_alphas_backward.csv'
            if molecula:
                his_mole = {'I_pos':[],'I_neg':[]}
                dfs = file[0]
                for i in range(len(file[0])):   
                    columns = dfs[i].columns
                    x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                    dx,dy = didv(x,y)
                    x_cut,y_cut = cut(x,y,xmin_hist ,xmax_hist  )
                   # I = round(integral(x_cut,y_cut),4)
                    x_cut_neg,y_cut_neg = cut(dx,dy,xmin,0-resolution )
                    I_neg = round(integral(x_cut_neg,y_cut_neg),4)
                    x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax_hist  )       
 
                    I_pos = round(integral(x_cut_pos,y_cut_pos),4)
                    his_mole['I_pos'].append(I_pos)
                    his_mole['I_neg'].append(I_neg)
                df = DataFrame(his_mole)
                df = df.to_csv(paste+'/'+'histogram_molecule_backward.csv')
                hist_path_molecula_backward = paste+'/'+'histogram_molecule_backward.csv'
            to_table(paste+'/',files_df,name=name_table)
            print("arquivos salvos na pasta "+ paste)

#Functions displayed for user

def Open_files():
    global folder
    global filenames
    global files_data
    global conca
    global stm
    stm = 'Nanosurf'
    gui = Tk()
    gui.geometry("600x500")
    gui.title('STS Data Analise for Omicron and Nanosurf STM by Rafael Reis Barreto')
    url = "https://i.postimg.cc/WzPSwVNh/logo.png"
    u = urlopen(url)
    raw_data = u.read()
    u.close()
    # In order to display the image in a GUI, you will use the 'PhotoImage' method of Tkinter. It will an image from the directory (specified path) and store the image in a variable.
    #icon = PhotoImage(file = "teste.png")
    icon = PhotoImage(data=raw_data)

    # Finally, to display the image you will make use of the 'Label' method and pass the 'image' variriable as a parameter and use the pack() method to display inside the GUI.
    label = Label(gui, image = icon)
    label.grid(row=0,column=1)


    label2 = Label(text='Which STM')
    label2.grid(row=1,column=1)

    def which():
        global stm
        stm = 'Nanosurf'
        print('Nanosurf mode selected')
        showinfo(
            title='Selected Mode',
            message='Nanosurf Mode Selected'
        )

    btnselec = ttk.Button(gui, text="Nanosurf Files", command=which)
    btnselec.grid(row=2,column=0)

    def which2():
        global stm
        stm = 'Omicron'
        print('Omicron mode selected')
        showinfo(
            title='Selected Mode',
            message='Omicron Mode Selected'
        )
    btnselec2 = ttk.Button(gui, text="Omicron Files", command=which2)
    btnselec2.grid(row=2,column=2)



    def getFolderPath():
                filetypes = (
            ('*.txt', '*.csv'),
            ('All files', '*.*')
        )

                global folder
                folder = filedialog.askopenfilename(filetypes=filetypes)
                print("File: ",folder)
                print("Uploaded")
                showinfo(
                title='Selected File',
                message=folder
            )
                
                gui.destroy()


    btnFind = ttk.Button(gui, text="Open File (Single File)", command=getFolderPath)
    btnFind.grid(row=3,column=0)

    def select_files():
        global filenames
        global conca
        filetypes = (
            ('*.txt', '*.csv'),
            ('All files', '*.*')
        )
        conca = False
        filenames = filedialog.askopenfilenames(
            title='Open files (Average)',
            #initialdir='/.',
            filetypes=filetypes)
        print("Files: ")
        for names in filenames:
            print(names)
        print("Uploaded")
        showinfo(
            title='Selected Files',
            message=filenames
        )
        gui.destroy()

    # open button
    open_button = ttk.Button(
        gui,
        text='Open files (Average)',
        command=select_files
    )

    #open_button.pack(expand=True)
    open_button.grid(row=3,column=1)

    def select_files2():
        global filenames
        global conca
        filetypes = (
            ('*.txt', '*.csv'),
            ('All files', '*.*')
        )
        conca = True
        filenames = filedialog.askopenfilenames(
            title='Open files (Concatenate)',
            #initialdir='/.',
            filetypes=filetypes)
        print("Files: ")
        for names in filenames:
            print(names)
        print("Uploaded")
        showinfo(
            title='Selected Files',
            message=filenames
        )
        gui.destroy()

    # open button
    open_button2 = ttk.Button(
        gui,
        text='Open files (Concatenate)',
        command=select_files2
    )

    #open_button.pack(expand=True)
    open_button2.grid(row=3,column=2)


    gui.mainloop()





    if filenames!=None:
        files_mult = []
        for item in filenames:
            files_mult.append(load_file(item,which=stm))
        
        for indice in range(len(filenames[0])):
            ct = len(filenames[0])-1-indice
            if filenames[0][ct]=='/':
                break
        folder = filenames[0][:ct+1] +'AAA'
        Vmin_glob =-9999;Vmax_glob =9999
        df_news = []
        if conca == False:
            for i in range(len(files_mult[0]['ixv'][0])):
                list_interps = []
                for j in range(len(files_mult)):
                    df = files_mult[j]['ixv'][0][i]
                    colum = df.columns
                    x_data = array(df[colum[0]]); y_data = array(df[colum[1]])
                    if x_data[0]>0:
                        xnew = [];ynew = []
                        for i in range(len(x_data)):
                            xnew.append(x_data[len(x_data)-i-1])
                            ynew.append(y_data[len(y_data)-1-i])
                        x_data = array(xnew)
                        y_data = array(ynew)

                    f = interpolate.interp1d(x_data,y_data)
                    list_interps.append(f)
                    if j ==0:
                        vmin = x_data.min();vmax = x_data.max()
                        if x_data.min()>Vmin_glob:
                            Vmin_glob = x_data.min()
                        if x_data.max()<Vmax_glob:
                            Vmax_glob = x_data.max()
                    else:
                        if x_data.min()>vmin:
                            vmin = x_data.min()
                        if x_data.max()<vmax:
                            vmax = x_data.max()
                        if x_data.min()>Vmin_glob:
                            Vmin_glob = x_data.min()
                        if x_data.max()<Vmax_glob:
                            Vmax_glob = x_data.max()
                V_new = arange(vmin,vmax,0.005)
                
                for k in range(len(list_interps)):
                    if k ==0:
                        ymedia = array(list_interps[k](V_new))
                    else:
                        ymedia=(ymedia+array(list_interps[k](V_new)))/2

                V_new = V_new
                df_media = DataFrame(V_new,columns = ['V'])
                df_media['I'+str(i)]=ymedia
                df_news.append(df_media)
                #printProgressBar(i,len(files_mult[0]['ixv'][0])-1)
######################
            if stm == 'Omicron':
                files_data = {'ixv':[df_news,['Average'],[Vmin_glob,Vmax_glob]],'ixv_idaevolta':[[],[],[]]}
            else:
                files_data = {'ixv':[df_news,files_mult[0]['ixv'][1],[Vmin_glob,Vmax_glob]],'ixv_idaevolta':[[],[],[]]}
            df_news = []
            if len(files_mult[0]['ixv_idaevolta'][0])!=0:
                for i in range(len(files_mult[0]['ixv_idaevolta'][0])):
                    list_interps = []
                    for j in range(len(files_mult)):
                        df = files_mult[j]['ixv_idaevolta'][0][i]
                        colum = df.columns
                        x_data = array(df[colum[0]]); y_data = array(df[colum[1]])
                        if x_data[0]>0:
                            xnew = [];ynew = []
                            for i in range(len(x_data)):
                                xnew.append(x_data[len(x_data)-i-1])
                                ynew.append(y_data[len(y_data)-1-i])
                            x_data = array(xnew)
                            y_data = array(ynew)

                        f = interpolate.interp1d(x_data,y_data)
                        list_interps.append(f)
                        if j ==0:
                            vmin = x_data.min();vmax = x_data.max()
                            if x_data.min()>Vmin_glob:
                                Vmin_glob = x_data.min()
                            if x_data.max()<Vmax_glob:
                                Vmax_glob = x_data.max()
                        else:
                            if x_data.min()>vmin:
                                vmin = x_data.min()
                            if x_data.max()<vmax:
                                vmax = x_data.max()
                            if x_data.min()>Vmin_glob:
                                Vmin_glob = x_data.min()
                            if x_data.max()<Vmax_glob:
                                Vmax_glob = x_data.max()
                    V_new = arange(vmin,vmax,0.005)
                    
                    for k in range(len(list_interps)):
                        if k ==0:
                            ymedia = array(list_interps[k](V_new))
                        else:
                            ymedia=(ymedia+array(list_interps[k](V_new)))/2

                    V_new = V_new
                    df_media = DataFrame(V_new,columns = ['V'])
                    df_media['I'+str(i)]=ymedia
                    df_news.append(df_media)
                    #printProgressBar(i,len(files_mult[0]['ixv'][0])-1)
                if stm == 'Omicron':
                    files_data['ixv_idaevolta'][0]=df_news
                    files_data['ixv_idaevolta'][1]=['Average_backward']
                    files_data['ixv_idaevolta'][2]=[Vmin_glob,Vmax_glob]
                else:
                    files_data['ixv_idaevolta'][0]=df_news
                    files_data['ixv_idaevolta'][1]=files_mult[0]['ixv'][1]
                    files_data['ixv_idaevolta'][2]=[Vmin_glob,Vmax_glob]

        elif conca == True:
            df_conca = []
            names_conca = []
            for j in range(len(files_mult)):
                if j ==0:
                    df = files_mult[j]['ixv'][0]
                    df_conca = df
                    if stm == 'Omicron':
                        names_conca =  files_mult[j]['ixv'][1]
                else:
                    df = files_mult[j]['ixv'][0]
                    df_conca=df_conca+df
                    if stm == 'Omicron':
                        names =  files_mult[j]['ixv'][1]
                        names_conca =names_conca +names
            for i in range(len(files_mult[0]['ixv'][0])):
                for j in range(len(files_mult)):
                    try:
                        df = files_mult[j]['ixv'][0][i]
                        colum = df.columns
                        x_data = array(df[colum[0]])
                        if x_data[0]>0:
                            x_data = sort(x_data)
                        if j ==0:
                            vmin = x_data.min();vmax = x_data.max()
                            if x_data.min()>Vmin_glob:
                                Vmin_glob = x_data.min()
                            if x_data.max()<Vmax_glob:
                                Vmax_glob = x_data.max()
                        else:
                            if x_data.min()>vmin:
                                vmin = x_data.min()
                            if x_data.max()<vmax:
                                vmax = x_data.max()
                            if x_data.min()>Vmin_glob:
                                Vmin_glob = x_data.min()
                            if x_data.max()<Vmax_glob:
                                Vmax_glob = x_data.max()
                    except IndexError:
                        pass
                #printProgressBar(i,len(files_mult[0]['ixv'][0])-1)
            print(len(df_conca))
            if stm == 'Omicron':
                files_data = {'ixv':[df_conca,names_conca,[Vmin_glob,Vmax_glob]],'ixv_idaevolta':[[],[],[]]}
            else:
                files_data = {'ixv':[df_conca,list(map(lambda x:str(x),range(len(df_conca)))),[Vmin_glob,Vmax_glob]],'ixv_idaevolta':[[],[],[]]}
    
    #############
            if len(files_mult[0]['ixv_idaevolta'][0])!=0:
                df_conca = []
                names_conca = []
                for j in range(len(files_mult)):
                    if j ==0:
                        df = files_mult[j]['ixv_idaevolta'][0]
                        df_conca = df
                        if stm == 'Omicron':
                            names_conca =  files_mult[j]['ixv_idaevolta'][1]
                    else:
                        df = files_mult[j]['ixv_idaevolta'][0]
                        df_conca=df_conca+df
                        if stm == 'Omicron':
                            names =  files_mult[j]['ixv_idaevolta'][1]
                            names_conca =names_conca +names
                for i in range(len(files_mult[0]['ixv_idaevolta'][0])):
                    for j in range(len(files_mult)):
                        try:
                            df = files_mult[j]['ixv_idaevolta'][0][i]
                            colum = df.columns
                            x_data = array(df[colum[0]])
                            if x_data[0]>0:
                                x_data = sort(x_data)
                            if j ==0:
                                vmin = x_data.min();vmax = x_data.max()
                                if x_data.min()>Vmin_glob:
                                    Vmin_glob = x_data.min()
                                if x_data.max()<Vmax_glob:
                                    Vmax_glob = x_data.max()
                            else:
                                if x_data.min()>vmin:
                                    vmin = x_data.min()
                                if x_data.max()<vmax:
                                    vmax = x_data.max()
                                if x_data.min()>Vmin_glob:
                                    Vmin_glob = x_data.min()
                                if x_data.max()<Vmax_glob:
                                    Vmax_glob = x_data.max()
                        except IndexError:
                            pass
                    #printProgressBar(i,len(files_mult[0]['ixv'][0])-1)
                print(len(df_conca))
                if stm == 'Omicron':
                    files_data['ixv_idaevolta'][0]=df_conca
                    files_data['ixv_idaevolta'][1]=names_conca
                    files_data['ixv_idaevolta'][2]=[Vmin_glob,Vmax_glob]
                else:
                    files_data['ixv_idaevolta'][0]=df_conca
                    files_data['ixv_idaevolta'][1]=list(map(lambda x:str(x),range(len(df_conca))))
                    files_data['ixv_idaevolta'][2]=[Vmin_glob,Vmax_glob]
    else:
        files_data = load_file(folder,which=stm)

def Display():
        #file = load_file(folder)

        if len(files_data['ixv_idaevolta'][0])!=0:
            file = files_data['ixv'];file_volta = files_data['ixv_idaevolta']
            def plot_curve(file,curve = 0,smooth = 5,delta = 10,resolution=0.010,volta= True,mm= False,xlim = (-1,1)):
                    global save_var
                    global curve_glob
                    global smooth_glob
                    global delta_glob
                    global selec 
                    curve_glob = curve; smooth_glob = smooth; delta_glob = delta
                    xmin,xmax = xlimSlider.value
                    curve = int(curve)
                    fig,ax= plt.subplots(1,2,figsize=(20,8))
                    dfs = file[0]
                    columns = dfs[curve].columns
                    x = dfs[curve][columns[0]];y = dfs[curve][columns[1]]*pow(10,9)

                    if mm:
                        x_cut,y_cut = cut(x,y,xmin,xmax )
                        a,b = regression(x_cut,y_cut)
                        a = round(a,3); b = round(b,3)
                        I = round(integral(x_cut,y_cut),3)
                        ax[0].plot(x_cut,x_cut*a+b,c='orange',label = 'a: ' + str(a) + ' b: '+str(b)+'\n'+'Integral: '+str(I))

                    p = int(smooth*len(y)/100)
                    if p%2==0:
                        p+=1
                        y = savgol_filter(y,p,1)
                    elif p==0:
                        pass
                    else:
                        y = savgol_filter(y,p,1)
                    if volta:
                        ax[0].plot(x,y,color = 'black',label = 'Forward curve')
                    else:   
                        ax[0].plot(x,y,color = 'black')
                        ax[0].legend()
                    if stm == 'Omicron':
                        ax[0].set_title('Curve from file: '+files_data['ixv'][1][curve])
                    ax[0].set_xlabel('Sample bias (V)')
                    ax[0].set_ylabel('Current (nA)')

                    dx,dy = didv(x,y)
                    if mm:
                        x_cut_neg,y_cut_neg = cut(dx,dy,xmin,0-resolution )
                        a_neg,b_neg = regression(x_cut_neg,y_cut_neg)
                        a_neg = round(a_neg,3); b_neg = round(b_neg,3)
                        I_neg = round(integral(x_cut_neg,y_cut_neg),3)
                        x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax )
                        a_pos,b_pos = regression(x_cut_pos,y_cut_pos)
                        a_pos = round(a_pos,3); b_pos = round(b_pos,3)                
                        I_pos = round(integral(x_cut_pos,y_cut_pos),3)
                        if volta:
                            ax[1].plot(dx,dy,color = 'black',label = 'F Int neg: '+str(I_neg)+' Int pos: '+str(I_pos))
                        else:
                            ax[1].plot(dx,dy,color = 'black',label = 'Int neg: '+str(I_neg)+' Int pos: '+str(I_pos))
                        ax[1].plot(x_cut_neg,a_neg*x_cut_neg+b_neg,c = 'orange',label = 'a_neg: ' + str(a_neg) + ' b_neg: '+str(b_neg))
                        ax[1].plot(x_cut_pos,a_pos*x_cut_pos+b_pos,c = 'gray',label = 'a_pos: ' + str(a_pos) + ' b_pos: '+str(b_pos))

                    gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                    if abs(round(typ,3))<=resolution:
                        tipo = 'neutro'
                    elif typ<-resolution:
                        tipo = 'n'
                    else:
                        tipo = 'p'
                    dyinterp = interpolate.interp1d(dx,dy)
                    ymin = dyinterp(xmin)
                    ymax = dyinterp(xmax)
                    ax[1].scatter([xmin,xmax],[ymin,ymax],s = 50, color = 'red')
                    if volta:
                        if mm == False:
                            ax[1].plot(dx,dy, label = 'F Gap '+ str(gap)+ ': Type ' + tipo,color = 'black')
                    else:
                        if mm == False:
                            ax[1].plot(dx,dy,color = 'black',label = 'Forward')
                    ax[1].set_title('Gap '+ str(gap)+ ': Type ' + tipo)
                    ax[1].set_xlabel('Sample bias (V)')
                    ax[1].set_ylabel('dI/dV (arb. units)')
                    ax[1].legend()    

                    if volta:
                        dfs = file_volta[0]
                        columns = dfs[curve].columns
                        x = dfs[curve][columns[0]];y = dfs[curve][columns[1]]*pow(10,9)
                        p = int(smooth*len(y)/100)
                        if p%2==0:
                            p+=1
                            y = savgol_filter(y,p,1)
                        elif p==0:
                            pass
                        else:
                            y = savgol_filter(y,p,1)

                        ax[0].plot(x,y,color ='red',label = 'Backward curve')
                        ax[0].set_xlabel('Sample bias (V)')
                        ax[0].set_ylabel('Current (nA)')
                        ax[0].legend()

                        dx,dy = didv(x,y)
                        gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                        if abs(round(typ,3))<=resolution:
                            tipo = 'neutro'
                        elif typ<-resolution:
                            tipo = 'n'
                        else:
                            tipo = 'p'
                        dyinterp = interpolate.interp1d(dx,dy)
                        ymin = dyinterp(xmin)
                        ymax = dyinterp(xmax)
                        ax[1].scatter([xmin,xmax],[ymin,ymax],s = 50, color = 'blue')
                        ax[1].plot(dx,dy, label = 'B Gap '+ str(gap)+ ': Type ' + tipo,color = 'red')
                        ax[1].set_xlabel('Sample bias (V)')
                        ax[1].set_ylabel('dI/dV (arb. units)')
                        ax[1].legend()    
            
            curve_slider=widgets.widgets.IntSlider(
            value=0,min=0,max=len(file[0])-1,	step=1,description='Select Curve: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            smoothingSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Smoothing: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            deltaSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            resolutionSlider = widgets.FloatText(
            value=0.010,
            description='Resolution of Doping Neutral (V):',
            disabled=False,layout=widgets.Layout(width='400px')	
    )
            voltaslider = widgets.Checkbox(value = True,description= 'Show the backward curve ',disable = False)

            molecu_metal_slider = widgets.Checkbox(value = False,description= 'Molecule or Metal ',disable = False)
            xlimSlider=widgets.FloatRangeSlider(value=[file[2][0], file[2][1]], min=file[2][0], max=file[2][1],  step=0.1, description='Range (V):',  disabled=False, continuous_update=False, orientation='horizontal',
            readout=True,  readout_format='.1f',layout=widgets.Layout(width='400px'))
            


            widgets.interact(plot_curve, file = widgets.fixed(file),  curve=curve_slider, smooth = smoothingSlider,delta = deltaSlider,resolution = resolutionSlider,
                             volta = voltaslider,mm = molecu_metal_slider,xlim = xlimSlider)
            #print(curve_glob,smooth_glob,delta_glob)


        else:
            file = files_data['ixv']
            def plot_curve(file,curve = 0,smooth = 5,delta = 10,resolution=0.010,mm=False,xlim = (-1,1)):
                    global save_var
                    global curve_glob
                    global smooth_glob
                    global delta_glob
                    global selec 
                    global stm
                    curve_glob = curve; smooth_glob = smooth; delta_glob = delta
                    xmin,xmax = xlimSlider.value
                    curve = int(curve)
                    fig,ax= plt.subplots(1,2,figsize=(20,8))
                    dfs = file[0]
                    columns = dfs[curve].columns
                    x = dfs[curve][columns[0]];y = dfs[curve][columns[1]]*pow(10,9)
                    if mm:
                        x_cut,y_cut = cut(x,y,xmin,xmax )
                        a,b = regression(x_cut,y_cut)
                        a = round(a,3); b = round(b,3)
                        I = round(integral(x_cut,y_cut),3)
                        ax[0].plot(x_cut,x_cut*a+b,c='orange',label = 'a: ' + str(a) + ' b: '+str(b)+'\n'+'Integral: '+str(I))

                    p = int(smooth*len(y)/100)
                    if p%2==0:
                        p+=1
                        y = savgol_filter(y,p,1)
                    elif p==0:
                        pass
                    else:
                        y = savgol_filter(y,p,1)

                    ax[0].plot(x,y,color = 'black')
                    ax[0].set_xlabel('Sample bias (V)')
                    ax[0].set_ylabel('Current (nA)')
                    if stm == 'Omicron':
                        ax[0].set_title('Curve from file: '+files_data['ixv'][1][curve])

                    dx,dy = didv(x,y)
                    if mm:
                        x_cut_neg,y_cut_neg = cut(dx,dy,xmin,0-resolution )
                        a_neg,b_neg = regression(x_cut_neg,y_cut_neg)
                        a_neg = round(a_neg,3); b_neg = round(b_neg,3)
                        I_neg = round(integral(x_cut_neg,y_cut_neg),3)
                        x_cut_pos,y_cut_pos = cut(dx,dy,0+resolution,xmax )
                        a_pos,b_pos = regression(x_cut_pos,y_cut_pos)
                        a_pos = round(a_pos,3); b_pos = round(b_pos,3)                
                        I_pos = round(integral(x_cut_pos,y_cut_pos),3)
                        ax[1].plot(dx,dy,color = 'black',label = 'Int neg: '+str(I_neg)+' Int pos: '+str(I_pos))
                        ax[1].plot(x_cut_neg,a_neg*x_cut_neg+b_neg,c = 'orange',label = 'a_neg: ' + str(a_neg) + ' b_neg: '+str(b_neg))
                        ax[1].plot(x_cut_pos,a_pos*x_cut_pos+b_pos,c = 'gray',label = 'a_pos: ' + str(a_pos) + ' b_pos: '+str(b_pos))
                        ax[1].legend()   


                    gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                    if abs(round(typ,3))<=resolution:
                        tipo = 'neutro'
                    elif typ<-resolution:
                        tipo = 'n'
                    else:
                        tipo = 'p'
                    dyinterp = interpolate.interp1d(dx,dy)
                    ymin = dyinterp(xmin)
                    ymax = dyinterp(xmax)
                    ax[1].scatter([xmin,xmax],[ymin,ymax],s = 50, color = 'red')
                    if mm == False:
                        ax[1].plot(dx,dy,color = 'black')
                    
                    ax[1].set_title('Gap '+ str(gap)+ ': Type ' + tipo)
                    ax[1].set_xlabel('Sample bias (V)')
                    ax[1].set_ylabel('dI/dV (arb. units)')
                     
        
            curve_slider=widgets.widgets.IntSlider(
            value=0,min=0,max=len(file[0])-1,	step=1,description='Select Curve: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            smoothingSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Smoothing: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            deltaSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold: ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
            resolutionSlider = widgets.FloatText(
            value=0.010,
            description='Resolution of Doping Neutral (V):',
            disabled=False,layout=widgets.Layout(width='400px')	
    )   
            molecu_metal_slider = widgets.Checkbox(value = False,description= 'Molecule or Metal ',disable = False)
            xlimSlider=widgets.FloatRangeSlider(value=[file[2][0], file[2][1]], min=file[2][0], max=file[2][1],  step=0.1, description='Range (V):',  disabled=False, continuous_update=False, orientation='horizontal',
            readout=True,  readout_format='.1f',layout=widgets.Layout(width='400px'))
            

            widgets.interact(plot_curve, file = widgets.fixed(file),  curve=curve_slider, smooth = smoothingSlider,delta = deltaSlider,resolution = resolutionSlider,mm = molecu_metal_slider,xlim = xlimSlider)
            #print(curve_glob,smooth_glob,delta_glob)

def Save_data():
    file = files_data['ixv']
    file_volta = files_data['ixv_idaevolta']
    for i in range(len(file[0])):
            selec[i] = True
    
    def str_to_int(A):
        A_new = []
        for item in A:
            try:
                A_new.append(int(item))
            except ValueError:
                pass
        return A_new

    def save_files(list_n,salvar = False,smooth = 5,delta = 5,resolution=0.01,save_um_file = False,metal= False, molecula = False,rang=(-1,1),volta = False):
        xmin,xmax = xlim.value
        global hist_path
        list_n = list_n.split(',')
        n_new = str_to_int(list_n)
        if salvar:
            print('Wainting ....')
            if volta:
                select_sts_ida(folder,file=file,n = n_new,smooth = smooth,delta = delta,resolution=resolution,
                       save_um_file = save_um_file,metal= metal, molecula = molecula,xmin_hist=  xmin,xmax_hist=  xmax,name_table ='forward' )
                select_sts_volta(folder,file=file_volta,n = n_new,smooth = smooth,delta = delta,resolution=resolution,
                       save_um_file = save_um_file,metal= metal, molecula = molecula,xmin_hist=  xmin,xmax_hist=  xmax,name_table ='backward')
                df_hist_ida = read_csv(hist_path_forwad)
                df_hist_volta = read_csv(hist_path_backward)
                df_hist = concat([df_hist_ida,df_hist_volta])
                df_hist = df_hist.reset_index(drop=True)
                hist_path= save_path+'/histogram.csv'
                df_hist.to_csv(hist_path)
                if metal:
                    df_hist_ida_metal = read_csv(hist_path_metal_forwad)
                    df_hist_volta_metal = read_csv(hist_path_metal_backward)
                    df_hist_metal = concat([df_hist_ida_metal,df_hist_volta_metal])
                    df_hist_metal = df_hist_metal.reset_index(drop=True)
                    hist_path_metal= save_path+'/histogram_alphas.csv'
                    df_hist_metal.to_csv(hist_path_metal)
                if molecula:
                    df_hist_ida_molecula = read_csv(hist_path_molecula_forwad)
                    df_hist_volta_molecula = read_csv(hist_path_molecula_backward)
                    df_hist_molecula = concat([df_hist_ida_molecula,df_hist_volta_molecula])
                    df_hist_molecula = df_hist_molecula.reset_index(drop=True)
                    hist_path_molecula= save_path+'/histogram_molecule.csv'
                    df_hist_molecula.to_csv(hist_path_molecula)      

            else:
                select_sts(folder,file=file,n = n_new,smooth = smooth,delta = delta,resolution=resolution,
                       save_um_file = save_um_file,metal= metal, molecula = molecula,xmin_hist=  xmin,xmax_hist=  xmax )
            print('Files saved')

    if len(files_data['ixv_idaevolta'][0])==0:   
        list_n = widgets.Textarea(
        value='',    placeholder="Write down the number of the files that you DON'T want to save. Separated by comma.",    
        description="List of curve numbers: ",    disabled=False, layout=widgets.Layout(width='400px'))
        print("Write down the number of the files that you DON'T want to save. Separated by comma.")
        save_buttom= widgets.ToggleButton(
        value=False,
        continuous_update=tuple,
        description='Save file and make a report',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Description',
        icon='check',
        layout=widgets.Layout(width='400px')
    )
        smoothingSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Smoothing (%): ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
        deltaSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold (%): ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
        resolutionSlider = widgets.FloatText(
        value=0.010,
        description='Resolution of Doping Neutral (V):',
        disabled=False,layout=widgets.Layout(width='400px')	
    )
        slidersave = widgets.Checkbox(value = False,description= 'Save files individually ',disable = False)

        metal = widgets.Checkbox(value = False,description= 'Save Histrogram metal ',disable = False)

        molecula = widgets.Checkbox(value = False,description= 'Save Histrogram molecule  ',disable = False)
        xlim=widgets.FloatRangeSlider(
        value=[file[2][0], file[2][1]],
        min=file[2][0],
        max=file[2][1],
        step=0.1,
        description='Range (V):',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',layout=widgets.Layout(width='400px')
    )

        widgets.interact(save_files, list_n = list_n,salvar = save_buttom,smooth = smoothingSlider,delta = deltaSlider,
                        resolution = resolutionSlider,save_um_file = slidersave,metal = metal,molecula = molecula,rang = xlim,volta = widgets.fixed(False))
    else:
        list_n = widgets.Textarea(
        value='',    placeholder="Write down the number of the files that you DON'T want to save. Separated by comma.",    
        description="List of curve numbers: ",    disabled=False, layout=widgets.Layout(width='400px'))
        print("Write down the number of the files that you DON'T want to save. Separated by comma.")
        save_buttom= widgets.ToggleButton(
        value=False,
        continuous_update=tuple,
        description='Save file and make a report',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Description',
        icon='check',
        layout=widgets.Layout(width='400px')
    )
        smoothingSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Smoothing (%): ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
        deltaSlider=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold (%): ',continuous_update=False,layout=widgets.Layout(width='400px')	)  
        resolutionSlider = widgets.FloatText(
        value=0.010,
        description='Resolution of Doping Neutral (V):',
        disabled=False,layout=widgets.Layout(width='400px')	
    )
        slidersave = widgets.Checkbox(value = False,description= 'Save files individually ',disable = False)

        metal = widgets.Checkbox(value = False,description= 'Save Histrogram metal ',disable = False)

        molecula = widgets.Checkbox(value = False,description= 'Save Histrogram molecule  ',disable = False)
        xlim=widgets.FloatRangeSlider(
        value=[file[2][0], file[2][1]],
        min=file[2][0],
        max=file[2][1],
        step=0.1,
        description='Range (V):',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',layout=widgets.Layout(width='400px')
    )
        volta = widgets.Checkbox(value = False,description= 'Save Backward files ',disable = False)


        widgets.interact(save_files, list_n = list_n,salvar = save_buttom,smooth = smoothingSlider,delta = deltaSlider,
                        resolution = resolutionSlider,save_um_file = slidersave,metal = metal,molecula = molecula,rang = xlim,volta = volta)
    
def Hist_plot():
    font = {'size'   : 14}
    matplotlib.rc('font', **font)
    if hist_path == 'not_path':
        gui = Tk()
        gui.geometry("400x400")
        gui.title('Nanosurf STS Data Analise by Rafael Reis Barreto')
        def getFolderPath():
            global hist_path
            hist_path = filedialog.askopenfilename()
            print("File: ",hist_path)
            print("Uploaded")
            gui.destroy()

        btnFind = ttk.Button(gui, text="Open a histogram file", command=getFolderPath)
        btnFind.grid(row=1,column=1)
        gui.mainloop()
    hist_file = read_csv(hist_path)


    def plot_hist(hist_file_col,bins,label,monocolor,rwidth):
            fig, ax = plt.subplots(figsize = (10,8))
            if monocolor:
                N,n_bins,patches =ax.hist(hist_file_col,bins = bins,rwidth= rwidth)
                # We'll color code by height, but you could use any scalar
                hist_file_col=array(hist_file_col)
                fracs = array(sorted(hist_file_col))/hist_file_col.max()

                # we need to normalize the data to 0..1 for the full range of the colormap
                norm = matplotlib.colors.Normalize(fracs.min(), fracs.max())

                # Now, we'll loop through our objects and set the color of each accordingly
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.viridis(norm(thisfrac))
                    thispatch.set_facecolor(color)
            else:
                ax.hist(hist_file_col,bins = bins,edgecolor='black',rwidth= rwidth)


            ax.set_xlabel(label)
            ax.set_ylabel('Counting')
            de = rwidth*(hist_file_col.max()-hist_file_col.min())/bins     
            print('Voltage width of each bar: %s (V)'%(round(de,3)))



    box_layout = widgets.Layout(
			border='dashed 1px gray',
			margin='0px 10px 10px 0px',
			padding='5px 5px 5px 5px',
			width='600px')


    style = {'description_width': 'initial'}

    panel=[{},{}]
    fig= []
    axes=[]

    
    for axis in [0,1]:
        panel[axis]['bins']=widgets.IntSlider(value=10,min=1,max=100,	
                                            step=1,description='Bins: ',continuous_update=False,layout=widgets.Layout(width='400px')	)
        panel[axis]['histcolor']=widgets.Checkbox(value = False,description= 'Hist color by x axis',disable = False)
        panel[axis]['rwidth']=widgets.FloatSlider(value=0.9,min=0.1,max=1,	
                                            step=.1,description='rwidth: ',continuous_update=False,layout=widgets.Layout(width='400px')	)

    colums = [hist_file['Gap(V)'],hist_file['Dop(Value)']]
    label = ['Gap (V)','Doping, Shifting from 0 (V) ']
    for axis in [0,1]:
        panel[axis]['output']=widgets.interactive_output(plot_hist,{'hist_file_col':widgets.fixed(colums[axis]),
                                                                    'bins':panel[axis]['bins'],
                                                                    'label': widgets.fixed(label[axis]),
                                                                    'monocolor': panel[axis]['histcolor'],
                                                                    'rwidth': panel[axis]['rwidth']})
        panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['bins'], panel[axis]['histcolor'],panel[axis]['rwidth']],layout=box_layout)
        panel[axis]['widget'].children[0].layout.height = '600px'


    outputPanel = widgets.HBox([panel[0]['widget'],panel[1]['widget']],layout=widgets.Layout(width='1200px'))
    return outputPanel
    
def Hist_metal():
    font = {'size'   : 14}
    matplotlib.rc('font', **font)
    if hist_path == 'not_path':
        gui = Tk()
        gui.geometry("400x400")
        gui.title('Nanosurf STS Data Analise by Rafael Reis Barreto')
        def getFolderPath():
            global hist_path_metal
            hist_path_metal = filedialog.askopenfilename()
            print("File: ",hist_path)
            print("Uploaded")
            gui.destroy()

        btnFind = ttk.Button(gui, text="Open a histogram file", command=getFolderPath)
        btnFind.grid(row=1,column=1)
        gui.mainloop()
    hist_file_metal = read_csv(hist_path_metal)
    def plot_hist(hist_file_col,bins,label,monocolor,rwidth):
            hist_file_col = array(hist_file_col)
            fig, ax = plt.subplots(figsize = (10,8))
            if monocolor:
                N,n_bins,patches =ax.hist(hist_file_col,bins = bins,rwidth= rwidth)
                # We'll color code by height, but you could use any scalar
                hist_file_col=array(hist_file_col)
                fracs = array(sorted(hist_file_col))/hist_file_col.max()

                # we need to normalize the data to 0..1 for the full range of the colormap
                norm = matplotlib.colors.Normalize(fracs.min(), fracs.max())

                # Now, we'll loop through our objects and set the color of each accordingly
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.viridis(norm(thisfrac))
                    thispatch.set_facecolor(color)
            else:
                ax.hist(hist_file_col,bins = bins,edgecolor='black',rwidth= rwidth)


            ax.set_xlabel(label)
            ax.set_ylabel('Counting')
            de = rwidth*(hist_file_col.max()-hist_file_col.min())/bins     
            print('Voltage width of each bar: %s (V)'%(round(de,3)))



    box_layout = widgets.Layout(
			border='dashed 1px gray',
			margin='0px 10px 10px 0px',
			padding='5px 5px 5px 5px',
			width='600px')


    style = {'description_width': 'initial'}

    panel=[{},{}]
    fig= []
    axes=[]

    
    for axis in [0,1]:
        panel[axis]['bins']=widgets.IntSlider(value=10,min=1,max=100,	
                                            step=1,description='Bins: ',continuous_update=False,layout=widgets.Layout(width='400px')	)
        panel[axis]['histcolor']=widgets.Checkbox(value = False,description= 'Hist color by x axis',disable = False)
        panel[axis]['rwidth']=widgets.FloatSlider(value=0.9,min=0.1,max=1,	
                                            step=.1,description='rwidth: ',continuous_update=False,layout=widgets.Layout(width='400px')	)

    colums = [hist_file_metal['a'],hist_file_metal['a_neg']+hist_file_metal['a_pos']]
    label = ['alpha IxV','alpha (dIxdV) ']
    for axis in [0,1]:
        panel[axis]['output']=widgets.interactive_output(plot_hist,{'hist_file_col':widgets.fixed(colums[axis]),
                                                                    'bins':panel[axis]['bins'],
                                                                    'label': widgets.fixed(label[axis]),
                                                                    'monocolor': panel[axis]['histcolor'],
                                                                    'rwidth': panel[axis]['rwidth']})
        panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['bins'], panel[axis]['histcolor'],panel[axis]['rwidth']],layout=box_layout)
        panel[axis]['widget'].children[0].layout.height = '600px'


    outputPanel = widgets.HBox([panel[0]['widget'],panel[1]['widget']],layout=widgets.Layout(width='1200px'))
    return outputPanel 

def Hist_molecule():
    font = {'size'   : 14}
    matplotlib.rc('font', **font)
    if hist_path_molecula == 'not_path':
        gui = Tk()
        gui.geometry("400x400")
        gui.title('Nanosurf STS Data Analise by Rafael Reis Barreto')
        def getFolderPath():
            global hist_path_molecula
            hist_path_molecula = filedialog.askopenfilename()
            print("File: ",hist_path)
            print("Uploaded")
            gui.destroy()

        btnFind = ttk.Button(gui, text="Open a histogram file", command=getFolderPath)
        btnFind.grid(row=1,column=1)
        gui.mainloop()
    hist_file_molecule = read_csv(hist_path_molecula)
    def plot_hist(hist_file_col,bins,label,rwidth,pizza):
            if pizza == True:
                fig, ax = plt.subplots(figsize = (10,8))
                coluna = hist_file_col.columns
                hist_file_cola = array(hist_file_col[coluna[0]])
                hist_file_colb = array(hist_file_col[coluna[1]])
                n_neg = 0;n_pos = 0;n_neuto=0
                for i in range(len(hist_file_cola)):
                    if abs(hist_file_cola[i])>abs(hist_file_colb[i]):
                        n_pos+=1
                    elif abs(hist_file_cola[i])<abs(hist_file_colb[i]):
                        n_neg+=1
                    else:
                        n_neuto +=0
                n = n_neg+n_neuto+n_pos
                labels = ['Int. Diff Neg.','Int. Diff. Zero','Int. Diff. Pos.']
                c = ['red','gray','blue']
                values = [n_neg*100/n,n_neuto*100,n_pos*100/n]
                wedges, texts, autote= ax.pie(values, autopct= '%1.1f%%',shadow=True,startangle=90,colors = c)
                ax.legend(wedges,labels, title='Int Values (%)',loc = 'center left')
            else:     
                coluna = hist_file_col.columns
                hist_file_cola = array(hist_file_col[coluna[0]])
                hist_file_colb = -array(hist_file_col[coluna[1]])
                fig, ax = plt.subplots(figsize = (10,8))
                ax.hist(hist_file_cola,bins = bins,edgecolor='black',rwidth= rwidth,color = 'blue')
                ax.hist(hist_file_colb,bins = bins,edgecolor='black',rwidth= rwidth,color = 'red')
                
                ax.set_xlabel(label)
                ax.set_ylabel('Counting')
                de = rwidth*(hist_file_cola.max()-hist_file_cola.min())/bins     
                print('Voltage width of each bar: %s (V)'%(round(de,3)))



    box_layout = widgets.Layout(
			border='dashed 1px gray',
			margin='0px 10px 10px 0px',
			padding='5px 5px 5px 5px',
			width='600px')


    style = {'description_width': 'initial'}

    panel=[{},{}]
    fig= []
    axes=[]

    
    for axis in [0,1]:
        panel[axis]['bins']=widgets.IntSlider(value=10,min=1,max=100,	
                                                step=1,description='Bins: ',continuous_update=False,layout=widgets.Layout(width='400px')	)

        panel[axis]['rwidth']=widgets.FloatSlider(value=0.9,min=0.1,max=1,	
                                                step=.1,description='rwidth: ',continuous_update=False,layout=widgets.Layout(width='400px')	)

    colums = [hist_file_molecule,hist_file_molecule]
    label = ['alpha IxV','Int. Values ']
    for axis in [0,1]:
        if axis ==0:
            panel[axis]['output']=widgets.interactive_output(plot_hist,{'hist_file_col':widgets.fixed(colums[axis]),
                                                                        'bins':widgets.fixed(panel[axis]['bins']),
                                                                        'label': widgets.fixed(label[axis]),
                                                                        'rwidth': widgets.fixed(panel[axis]['rwidth']),
                                                                        'pizza':widgets.fixed(True),
                                                                        })
            panel[axis]['widget'] = widgets.VBox([panel[axis]['output']],layout=box_layout)
            panel[axis]['widget'].children[0].layout.height = '600px'      
        else:
            panel[axis]['output']=widgets.interactive_output(plot_hist,{'hist_file_col':widgets.fixed(colums[axis]),
                                                                        'bins':panel[axis]['bins'],
                                                                        'label': widgets.fixed(label[axis]),
                                                                        'rwidth': panel[axis]['rwidth'],
                                                                        'pizza':widgets.fixed(False),
                                                                        })
            panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['bins'],panel[axis]['rwidth']],layout=box_layout)
            panel[axis]['widget'].children[0].layout.height = '600px'


    outputPanel = widgets.HBox([panel[0]['widget'],panel[1]['widget']],layout=widgets.Layout(width='1200px'))
    return outputPanel 


def Grid_plot():
    font = {'size'   : 14}
    matplotlib.rc('font', **font)
    #file = load_file(folder)
    def save_matrix(x,y,M,name,tipo):
            folder_save = 'sts_saves'
            for i in range(1,len(name[:-4])):
                if name[len(name[:-4]) -i] == '/':
                    ct = len(name[:-4]) -i
                    break
            try:
                #folder_name = path.join(folder_save,path_file[:-4][ct:])
                folder_name= folder_save+name[:-4][ct:]
            except UnboundLocalError:
                ct = 0
                #folder_name = path.join(folder_save,path_file[ct:-4])
                folder_name= folder_save+'/'+ name[:-4][ct:]
            try: 
                mkdir(folder_save )
            except FileExistsError:
                pass
                
            paste =folder_name
            try:
                mkdir(paste)
            except FileExistsError:
                pass
            arq = open(paste+'/grid'+tipo+'_.txt','w')
            arq.writelines('Grid of '+tipo)
            line ='x: '
            for item in x:
                line+=str(item)+','
            line +='\n'
            arq.writelines(line)
            line ='y: '
            for item in y:
                line+=str(item)+','
            line +='\n'            
            arq.writelines(line)
            line = 'z :'
            for i in range(len(M)):
                if i ==0:
                    for j in range(len(M[0])):
                        line+=str(M[i][j])+','
                    line += '\n'
                    arq.writelines(line)
                else:
                    line =''
                    for j in range(len(M[0])):
                        line+=str(M[i][j])+','
                    line += '\n'
                    arq.writelines(line)
            arq.close()
            print('File saved at '+paste)

    if stm =='Nanosurf':
        file = files_data['ixv']
        def get_matrix(dfs,delta,nx=20,V=0):
            files = []
            for i in range(len(dfs)):
                columns = dfs[i].columns
                x = dfs[i][columns[0]];y = dfs[i][columns[1]]*pow(10,9)
                files.append([x,y])
            ny = int(len(dfs)/nx)
            M_didv = [];M_dop = [];M_gap = []
            for i in range(nx):
                col_didv = [];col_dop = [];col_gap = []
                for j in range(ny):
                    x,y = files[j+i*ny]
                    dx,dy = didv(x,y)
                    gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                    dyinterp = interpolate.interp1d(dx,dy)

                    col_didv.append(dyinterp(V))
                    col_dop.append(round(typ,3))
                    col_gap.append(round(gap,3))
                M_didv.append(array(col_didv));M_dop.append(array(col_dop));M_gap.append(array(col_gap))

            M_didv = array(M_didv);M_dop = array(M_dop);M_gap = array(M_gap)
            x = arange(nx);y=arange(ny)
            return [x,y,[M_didv,M_gap,M_dop]]

        def plot_mapa2(file,whichmap,delta = 5,interpolation= False,mult=2,cmap = 'viridis',nx = 20,save_mat=False):
            x,y,mapas = get_matrix(file,delta,nx=nx)
            M=mapas[whichmap]

            fig,ax= plt.subplots(figsize=(10,8))

            if interpolation == True:
                nx = len(x);ny=len(y)

                x,y,M = inter_x(x,y,M,dx=mult*nx,dy = mult*ny)

            Y,X = meshgrid(y,x)
            im = ax.pcolormesh(X,Y,M,cmap = cmap)
            if save_mat:
                if axis==1:
                    save_matrix(x,y,M,folder,' Gap')
                elif axis==2:
                    save_matrix(x,y,M,folder,' Doping')
            if whichmap==1:
                cbar = fig.colorbar(im)
                cbar.set_ticks([M.min(),M.max()],labels= ['Gap ','No Gap'])
            elif whichmap==2:
                cbar = fig.colorbar(im)
                cbar.set_ticks([M.min(),0,M.max()],labels= ['N ','Neutro','P'])

            #im2 = ax[2].pcolormesh(X,Y,M_dop,cmap = cmap)
            #cbar = fig.colorbar(im2)
            #cbar.set_ticks([M_dop.min(),0,M_dop.max()],labels= ['N ','Neutro','P'])

        def plot_mapa(file,whichmap,V = 0,delta = 5,interpolation= False,mult=2,cmap = 'viridis',nx =20,save_mat=False):
            x,y,mapas = get_matrix(file,delta,V=V,nx = nx)
            M=mapas[whichmap]

            fig,ax= plt.subplots(figsize=(10,8))

            if interpolation == True:
                nx = len(x);ny=len(y)

                x,y,M = inter_x(x,y,M,dx=mult*nx,dy = mult*ny)
            if save_mat:
                save_matrix(x,y,M,folder,' dIdV at '+str(V)+' V')
            Y,X = meshgrid(y,x)
            im = ax.pcolormesh(X,Y,M,cmap = cmap)
            cbar = fig.colorbar(im)


        box_layout = widgets.Layout(
                border='dashed 1px gray',
                margin='0px 10px 10px 0px',
                padding='5px 5px 5px 5px',
                width='600px')

        panel=[{},{},{}]
        fig= []
        axes=[]
        for axis in [0,1,2]:
            if axis ==0:
                panel[axis]['vslider']=widgets.widgets.FloatSlider(
                value=0,min=file[2][0],max=file[2][1],	step=.2,description='Sample Bias (V): ',continuous_update=False,layout=widgets.Layout(width='300px')	)
            panel[axis]['deltaslider']=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold: ',continuous_update=False,layout=widgets.Layout(width='300px')	)  
            panel[axis]['interpolation']=widgets.Checkbox(value = False,description= 'Interpolation',disable = False)
            panel[axis]['mult']=widgets.widgets.IntSlider(
            value=2,min=1,max=20,	step=1,description='Multiplicity of interpolation ',continuous_update=False,layout=widgets.Layout(width='300px'))	
            panel[axis]['nx']=widgets.IntText(  value=20, description='Number of points per line. Nanosurf software',  disabled=False
    )
            panel[axis]['colormap'] = widgets.Dropdown(
            options=['bone_r', 'inferno', 'viridis','plasma', 'cividis','gray','OrRd','PuBuGn','coolwarm','bwr','terrain'],
            value='viridis',
            description='Colormap:',
            )
            panel[axis]['save']  = widgets.Checkbox(value = False,description= 'Save matrix file ',disable = False)

        for axis in [0,1,2]:

            if axis ==0:
                panel[axis]['output']=widgets.interactive_output(plot_mapa,{'file':widgets.fixed(file[0]),
                                                                        'whichmap':widgets.fixed(axis),
                                                                        'V':panel[axis]['vslider'],
                                                                        'delta': panel[axis]['deltaslider'],
                                                                        'interpolation': panel[axis]['interpolation'],
                                                                        'mult': panel[axis]['mult'],
                                                                        'cmap':panel[axis]['colormap'],
                                                                        'nx':panel[axis]['nx'],'save_mat':panel[axis]['save']
                                                                        })
                panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['vslider'], panel[axis]['deltaslider'],
                                                panel[axis]['interpolation'],panel[axis]['mult'],panel[axis]['colormap'],panel[axis]['nx'],panel[axis]['save'] ],layout=box_layout)
                panel[axis]['widget'].children[0].layout.height = '400px'
            else:
                panel[axis]['output']=widgets.interactive_output(plot_mapa2,{'file':widgets.fixed(file[0]),
                                                                        'whichmap':widgets.fixed(axis),
                                                                        'delta': panel[axis]['deltaslider'],
                                                                        'interpolation': panel[axis]['interpolation'],
                                                                        'mult': panel[axis]['mult'],
                                                                        'cmap':panel[axis]['colormap'],
                                                                        'nx':panel[axis]['nx'],
                                                                        'save_mat':panel[axis]['save']
                                                                        })
                panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['deltaslider'],
                                                panel[axis]['interpolation'],panel[axis]['mult'],panel[axis]['colormap'],panel[axis]['nx'],panel[axis]['save'] ],layout=box_layout)
            panel[axis]['widget'].children[0].layout.height = '600px'


        outputPanel = widgets.HBox([panel[0]['widget'],panel[1]['widget'],panel[2]['widget']],layout=widgets.Layout(width='1800px'))
        return outputPanel
    elif stm == 'Omicron':

        def map_spip_temp(map_file):
            arquivo = open_file(map_file)
            for line in arquivo:
                if ('z-range' in line) == True:
                    z_range = float(line[3])*pow(10,-9)
                elif ('z-points' in line) == True:
                    z_points = int(line[3])                 
                

            V_start = -z_range/2
            V_end = z_range/2
            dv = (V_end -V_start)/(z_points)
            V = []
            number = V_start
            for i in range(z_points):
                V.append(number)
                number+=dv
                number  = round(number,3)
            V = array(V)

            return [arquivo,V]


        def map_spip(arquivo,delta = 5,V_thre=0):
            matriz = []
            i = 0
            for line in arquivo:
                if ('x-pixels' in line) == True:
                    n = int(line[3])
                elif ('y-pixels' in line) == True:
                    m = int(line[3])
                elif ('x-length' in line) == True:
                    x_leng = float(line[3]) 
                elif ('y-length' in line) == True:
                    y_leng = float(line[3])
                elif ('x-offset' in line) == True:
                    x_ofsset = float(line[3]) 
                elif ('y-offset' in line) == True:
                    y_ofsset = float(line[3]) 
                elif ('z-range' in line) == True:
                    z_range = float(line[3])*pow(10,-9)
                elif ('z-points' in line) == True:
                    z_points = int(line[3])                 
                
                elif  ('Start' in line) == True and ('of' in line) == True and ('Data:' in line) == True:
                    break 
                i+=1

            
            dx = (x_leng)/n
            x = arange(-x_leng/2,x_leng/2,dx)-x_ofsset
            dy = (y_leng)/m
            y = arange(-y_leng/2,y_leng/2,dy)-y_ofsset


        
            sts_all =array(list(map(lambda x:array(list(map( lambda y:float(y),x))), arquivo[i+1:])))
            V_points = len(sts_all[0])
            V_start = -z_range/2
            V_end = z_range/2
            dv = (V_end -V_start)/(V_points)
            V = []
            number = V_start
            for i in range(V_points):
                V.append(number)
                number+=dv
                number  = round(number,3)
            V = array(V)
            matriz = []
            for j in range(m):
                col = []
                for i in range(n):
                    col.append(sts_all[j+(n-i-1)*n])
                matriz.append(array(col))

            nx = len(x);ny =len(y)
            M_didv = [];M_dop=[];M_gap = []
            for i in range(nx):
                col_didv = [];col_dop = [];col_gap = []
                for j in range(ny):

                    dx,dy = didv(V,matriz[i][j])
                    gap,typ,xmin,xmax = gap_type(dx,dy,delta/100)
                    
                    dyinterp = interpolate.interp1d(dx,dy)

                    col_didv.append(dyinterp(V_thre))
                    col_dop.append(round(typ,3))
                    col_gap.append(round(gap,3))
                M_didv.append(array(col_didv));M_dop.append(array(col_dop));M_gap.append(array(col_gap))

            M_didv = array(M_didv);M_dop = array(M_dop);M_gap = array(M_gap)
            return [x,y,[M_didv,M_gap,M_dop]]

        
        def plot_mapa2(arquivo,whichmap,delta = 5,interpolation= False,mult=2,cmap = 'viridis',save_mat=False):
            global folder
            x,y,mapas =map_spip(arquivo,delta = delta)
            M=mapas[whichmap]

                
            fig,ax= plt.subplots(figsize=(10,8))

            if interpolation == True:
                nx = len(x);ny=len(y)

                x,y,M = inter_x(x,y,M,dx=mult*nx,dy = mult*ny)

            Y,X = meshgrid(y,x)
            im = ax.pcolormesh(X,Y,M,cmap = cmap)
            ax.set_xlabel('Distance (nm)')
            ax.set_ylabel('Distance (nm)')
            if whichmap==1:
                cbar = fig.colorbar(im)
                cbar.set_ticks([M.min(),M.max()],labels= ['Gap ','No Gap'])
            elif whichmap==2:
                cbar = fig.colorbar(im)
                cbar.set_ticks([M.min(),0,M.max()],labels= ['N ','Neutro','P'])

            if save_mat:
                if axis==1:
                    save_matrix(x,y,M,folder,' Gap')
                elif axis==2:
                    save_matrix(x,y,M,folder,' Doping')
            #im2 = ax[2].pcolormesh(X,Y,M_dop,cmap = cmap)
            #cbar = fig.colorbar(im2)
            #cbar.set_ticks([M_dop.min(),0,M_dop.max()],labels= ['N ','Neutro','P'])

        def plot_mapa(arquivo,whichmap,V = 0,delta = 5,interpolation= False,mult=2,cmap = 'viridis',save_mat=False):
            x,y,mapas =map_spip(arquivo,delta = delta,V_thre=V)
            M=mapas[whichmap]
            fig,ax= plt.subplots(figsize=(10,8))

            if interpolation == True:
                nx = len(x);ny=len(y)

                x,y,M = inter_x(x,y,M,dx=mult*nx,dy = mult*ny)

            Y,X = meshgrid(y,x)
            im = ax.pcolormesh(X,Y,M,cmap = cmap)
            cbar = fig.colorbar(im)
            ax.set_xlabel('Distance (nm)')
            ax.set_ylabel('Distance (nm)')
            if save_mat:
                save_matrix(x,y,M,folder,' dIdV at '+str(V)+' V')



        box_layout = widgets.Layout(
                border='dashed 1px gray',
                margin='0px 10px 10px 0px',
                padding='5px 5px 5px 5px',
                width='600px')

        arquivo,V = map_spip_temp(folder)
        
        
        panel=[{},{},{}]
        fig= []
        axes=[]
        for axis in [0,1,2]:
            if axis ==0:
                panel[axis]['vslider']=widgets.widgets.FloatSlider(
                value=0,min=V.min(),max=V.max(),	step=.2,description='Sample Bias (V): ',continuous_update=False,layout=widgets.Layout(width='300px')	)
            panel[axis]['deltaslider']=widgets.widgets.FloatSlider(
            value=5,min=.5,max=20,	step=.5,description='Threshold: ',continuous_update=False,layout=widgets.Layout(width='300px')	)  
            panel[axis]['interpolation']=widgets.Checkbox(value = False,description= 'Interpolation',disable = False)
            panel[axis]['mult']=widgets.widgets.IntSlider(
            value=2,min=1,max=20,	step=1,description='Multiplicity of interpolation ',continuous_update=False,layout=widgets.Layout(width='300px'))	
            panel[axis]['colormap'] = widgets.Dropdown(
            options=['bone_r', 'inferno', 'viridis','plasma', 'cividis','gray','OrRd','PuBuGn','coolwarm','bwr','terrain'],
            value='viridis',
            description='Colormap:',
            )       
            panel[axis]['save']  = widgets.Checkbox(value = False,description= 'Save matrix file ',disable = False)


        for axis in [0,1,2]:

            if axis ==0:
                panel[axis]['output']=widgets.interactive_output(plot_mapa,{'arquivo':widgets.fixed(arquivo),
                                                                        'whichmap':widgets.fixed(axis),
                                                                        'V':panel[axis]['vslider'],
                                                                        'delta': panel[axis]['deltaslider'],
                                                                        'interpolation': panel[axis]['interpolation'],
                                                                        'mult': panel[axis]['mult'],
                                                                        'cmap':panel[axis]['colormap'],
                                                                        'save_mat':panel[axis]['save']
                                                                        })
                panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['vslider'], panel[axis]['deltaslider'],
                                                panel[axis]['interpolation'],panel[axis]['mult'],panel[axis]['colormap'],panel[axis]['save']],layout=box_layout)
                panel[axis]['widget'].children[0].layout.height = '400px'
            else:
                panel[axis]['output']=widgets.interactive_output(plot_mapa2,{'arquivo':widgets.fixed(arquivo),
                                                                        'whichmap':widgets.fixed(axis),
                                                                        'delta': panel[axis]['deltaslider'],
                                                                        'interpolation': panel[axis]['interpolation'],
                                                                        'mult': panel[axis]['mult'],
                                                                        'cmap':panel[axis]['colormap'],
                                                                        'save_mat':panel[axis]['save']
                                                                        })
                panel[axis]['widget'] = widgets.VBox([panel[axis]['output'],panel[axis]['deltaslider'],
                                                panel[axis]['interpolation'],panel[axis]['mult'],panel[axis]['colormap'],panel[axis]['save']],layout=box_layout)
            panel[axis]['widget'].children[0].layout.height = '600px'


        outputPanel = widgets.HBox([panel[0]['widget'],panel[1]['widget'],panel[2]['widget']],layout=widgets.Layout(width='1800px'))
        return outputPanel




