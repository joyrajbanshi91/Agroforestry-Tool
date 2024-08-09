from tkinter import *
from tkinter.filedialog import askopenfilename, asksaveasfilename
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pylab import rcParams
#rcParams['figure.figsize'] = 10,5

###### Sigmoid difference to calculate the carbon sequestration for initial forest area

def btn_clicked():
    try:
        float(entry0.get())
        float(entry1.get())
        float(entry2.get())
        float(entry3.get())
        float(entry4.get())
        float(entry5.get())
    except ValueError:
        label.config(text="You have entered an invalid input")
    
        
        
    Ini_F_Veg_C=eval(entry0.get())
    Ini_F_M_Age=eval(entry1.get())
    Ini_F_M_Age_1=eval(entry2.get())
    Ini_F_Soil_C=eval(entry3.get())
    soil_Mature_Age=eval(entry4.get())
    Rotational_Year=eval(entry5.get())
    
    
    def sig_diff_(start_year,end_year,maturity):
        start_year=start_year
        end_year=end_year
        Mature_Age=maturity
        sig_diff=[]
        prv_sig=( 1 - math.exp( ( -3.0 * 0 ) / Mature_Age ))**2
        for i in range(start_year,end_year+1,1):
            offyr=i-start_year
            curr_sig=(1 - math.exp( ( -3.0 * (offyr+1 ) ) / Mature_Age ))**2
            sig_diff.append(curr_sig-prv_sig)
            prv_sig=curr_sig
        sig_diff=np.array(sig_diff)
        return sig_diff
#### Sigmoid difference to calculate the carbon sequestration for the additional forest area
    def sig_diff(start_year,end_year,carbon_diff,maturity):
        start_year=start_year
        end_year=end_year
        carbon_diff=carbon_diff
        Mature_Age=maturity
        sig_diff=[]
        prv_sig=( 1 - math.exp( ( -3.0 * 0 ) / Mature_Age ))**2
        for i in range(start_year,end_year+1,1):
            offyr=i-start_year
            curr_sig=(1 - math.exp( ( -3.0 * (offyr+1 ) ) / Mature_Age ))**2
            sig_diff.append(curr_sig-prv_sig)
            prv_sig=curr_sig
        sig_diff=np.array(sig_diff)
        return sig_diff

########## soil Emission function for the initial forest area 
    def Soil_exponential_Init(start_year,end_year,timescale):
        start_year=start_year 
        end_year=end_year
        SoilTimeScale=timescale  
        halfLife=SoilTimeScale/10
        log2 = math.log( 2.0)
        Constant_K = log2 / halfLife;      
        cumStockDiff_t1 = 0
        yearCounter = 0
        exponential_diff=[]
        for i in range(start_year,end_year+1,1):
            yearCounter+=1
            cumStockDiff_t2 = ( 1.0 - math.exp( -1.0 * Constant_K * yearCounter ))
            exponential_diff.append(cumStockDiff_t2-cumStockDiff_t1)
            cumStockDiff_t1 = cumStockDiff_t2
        exponential_diff=np.array(exponential_diff)
        return exponential_diff

########## soil Emission function for the additional forest area 

    def Soil_exponential(start_year,end_year,carbon_diff,timescale):
        start_year=start_year 
        end_year=end_year
        carbon_diff=carbon_diff
        SoilTimeScale=timescale  
        halfLife=SoilTimeScale/10
        log2 = math.log( 2.0)
        Constant_K = log2 / halfLife;      
        cumStockDiff_t1 = 0
        yearCounter = 0
        exponential_diff=[]
        for i in range(start_year,end_year+1,1):
            yearCounter+=1
            cumStockDiff_t2 = ( 1.0 - math.exp( -1.0 * Constant_K * yearCounter ))
            exponential_diff.append(cumStockDiff_t2-cumStockDiff_t1)
            cumStockDiff_t1 = cumStockDiff_t2
        exponential_diff=np.array(exponential_diff)*carbon_diff
        return exponential_diff
#########################################  Main Code Starts from Here... ##############
    label.config(text="Please Select the input file")    
    Openfile=askopenfilename(filetypes=[("csv file","*.csv")],defaultextension=".csv")
    df = pd.read_csv(Openfile)
    Veg_C , Soil_C , Mature_Age = df['Veg_c'].iloc[0] , df['Soil_c'].iloc[0] , df['Mature_age'].iloc[0]
    df.drop(['Veg_c', 'Soil_c','Mature_age'], axis=1, inplace=True)

    Year = np.arange(df["Year"].iloc[0],df["Year"].iloc[-1]+1,1)
    df2 = pd.DataFrame(Year, columns=['Year'])

    merged_data= df2.merge(df, how='outer',on=["Year"])
    merged_data['LU_Area_ha']=merged_data['LU_Area_ha'].interpolate(method='linear')

    merged_data['land_diff'] = merged_data['LU_Area_ha'].diff(periods=-1)
    merged_data['land_diff']= merged_data['land_diff'].shift(periods=1, fill_value=0)

    merged_data['Abv_g_cd']=merged_data['land_diff']*Veg_C
    merged_data['Below_g_cd']=merged_data['land_diff']*Soil_C
    
    ######################## Calculation Existing Forest Carbon Sequestration
        
    sig_diff_val=sig_diff(merged_data['Year'][1],merged_data['Year'][1]+1000,1,Ini_F_M_Age) # Here 1000 is the years from the starting years
    Ini_forest=pd.DataFrame(sig_diff_val,columns=['sig_diff'])    
    Ini_F_Area=merged_data['LU_Area_ha'].iloc[0]
    Ini_F_C_density=Ini_F_Area*Ini_F_Veg_C
    Ini_forest['Veg_C_Seq']=Ini_forest['sig_diff']*Ini_F_C_density
    Ini_forest_C_Seq=Ini_forest['Veg_C_Seq'].iloc[Ini_F_M_Age_1-1 : Ini_F_M_Age_1+(len(merged_data))-1]
    Ini_forest_C_Seq=np.array(-abs(Ini_forest_C_Seq))

    
    ################################ Calculation of Existing Soil Carbon Sequestration #########################

    expo_diff = Soil_exponential_Init(merged_data['Year'][1],merged_data['Year'][1]+1000,soil_Mature_Age) # Here 1000 is the years from the starting years
    expo_soil_forest = pd.DataFrame(expo_diff, columns=['expo_diff'])
    Ini_F_Area=merged_data['LU_Area_ha'].iloc[0]
    Ini_F_Soil_C_density=Ini_F_Area*Ini_F_Soil_C
    expo_soil_forest['Soil_C_Seq']=expo_soil_forest['expo_diff']*Ini_F_Soil_C_density
    Ini_forest_Soil_C_Seq=expo_soil_forest['Soil_C_Seq'].iloc[soil_Mature_Age-1 : soil_Mature_Age+(len(merged_data))-1]
    Ini_forest_Soil_C_Seq=np.array(-abs(Ini_forest_Soil_C_Seq)) 
    
    ###########################   Calculation of Additional Vegetation Carbon Emissions ###############################

    vec=[]
    vec2=[]
    c=0
    for i,j in zip (merged_data['Abv_g_cd'],merged_data['Below_g_cd']):
        sy=merged_data['Year'][c]
        endy=merged_data['Year'][len(merged_data)-1]
        if i==0:
            vec.append (i)
        elif i<0:
            a=(sig_diff_(sy,endy,Mature_Age))
            b=(Soil_exponential(sy, endy,j,soil_Mature_Age))
            vec.append((a))
            vec2.append(b)
        elif i>0:
            a=(sig_diff_(sy,endy,Mature_Age))
            b=(Soil_exponential(sy, endy,j,soil_Mature_Age))[0]
            vec.append(a)
            vec2.append(b)
        c+=1
    ###########################  Calculation of Additional  Soil carbon Emissions #################
    first_year=merged_data['Year'][1]
    Soil_df=pd.DataFrame({'Year_{}'.format(first_year):vec2[0]})
    c=1
    for i in (merged_data['Year'].iloc[2 :]):
        Soil_df.loc[:, 'Year_{}'.format(i)] =pd.Series(vec2[c])
        c+=1

    c=0
    i = merged_data['Year'][1]
    for j in Soil_df.columns:
        Soil_df["Year_{}".format(i)] = Soil_df[j].shift(periods=c, fill_value=0)
        i += 1
        c+=1
    Soil_df=Soil_df.fillna(0)
    Soil_df["Soil_Emission"]=Soil_df.sum(axis=1)
        

    ################## Creating the final Vegetation Emissions  ############################################

    first_year=merged_data['Year'][1]
    new_df=pd.DataFrame({'Year_{}'.format(first_year):vec[1]})
    c=2
    for i in (merged_data['Year'].iloc[2 :]):
        new_df.loc[:, 'Year_{}'.format(i)] =pd.Series(vec[c])
        c+=1

    c=0
    i = merged_data['Year'][1]
    for j in new_df.columns:
        new_df["Year_{}".format(i)] = new_df[j].shift(periods=c, fill_value=0)
        i += 1
        c+=1
    new_df=new_df.fillna(0)

    Prev_Density=Veg_C
    Prev_Stock = merged_data['LU_Area_ha'].iloc[0]*Ini_F_Veg_C
    merged_data.drop(index=merged_data.index[0],axis=0,inplace=True)
    Abv_c_df = pd.DataFrame(columns=['Year'])
    Abv_G_Stock=[]
    c=1
    Total_Abv_c=[]
    for i in range(len(merged_data)):
        carbon_diff=merged_data['land_diff'].iloc[i]*Prev_Density

        if carbon_diff<0: 
            carbon_diff=merged_data['land_diff'].iloc[i]*Veg_C 
            C_seq_rate=carbon_diff*new_df.iloc[:,i]
            Abv_c_df[new_df.columns[i]]=C_seq_rate
            Abv_c_df=Abv_c_df.fillna(0)
            C_seq=Abv_c_df.sum(axis=1)[i]
            Total_Abv_c_Seq=Ini_forest_C_Seq[i]+C_seq
            Total_Abv_c.append(Total_Abv_c_Seq)
            Current_Stock=Prev_Stock+abs(Total_Abv_c_Seq)
            Current_C_density=Current_Stock/merged_data['LU_Area_ha'][c]
            Abv_G_Stock.append(Current_Stock)
            
        elif carbon_diff>0:   
            futureadj= merged_data['land_diff'][c]*(Veg_C-(Prev_Density))
            C_seq_rate=futureadj*new_df.iloc[:,i]
            ## Possible change should be done here to correct future adjustment 
            Abv_c_df[new_df.columns[i]]=C_seq_rate
            Abv_c_df=Abv_c_df.fillna(0)
            C_seq=Abv_c_df.sum(axis=1)[i]
            Total_Abv_c_Seq=carbon_diff+Ini_forest_C_Seq[i]+C_seq
            Total_Abv_c.append(Total_Abv_c_Seq)
            Current_Stock=Prev_Stock-(Total_Abv_c_Seq)
            Current_C_density=Current_Stock/merged_data['LU_Area_ha'][c]
            Abv_G_Stock.append(Current_Stock)
            
        Prev_Stock=Current_Stock
        Prev_Density=Current_C_density        
        c+=1
    ############################  Calculation of Agroforstry Rotational carbon Emissions #################
    Abv_c_df = Abv_c_df.drop('Year', axis=1)
    

    In_rot_yr=Rotational_Year
    if In_rot_yr>0:
        c=0
        rep=[]
        for col in Abv_c_df.columns:
            if In_rot_yr!= len(Abv_c_df):
                Repeat_array=np.array(Abv_c_df[col])[c:In_rot_yr]
                No_Repeat=int(len(Abv_c_df))
                Repeat_array_f=np.tile(Repeat_array,No_Repeat)
                Repeat_array_f=Repeat_array_f[0:len(Abv_c_df)]
                rep.append(Repeat_array_f)
            elif In_rot_yr>= len(Abv_c_df):
                Repeat_array=np.array(Abv_c_df[col])[c:len(Abv_c_df)]
                No_Repeat=int(len(Abv_c_df))
                Repeat_array_f=np.tile(Repeat_array,No_Repeat)
                Repeat_array_f=Repeat_array_f[0:len(Abv_c_df)]
                rep.append(Repeat_array_f)
            c+=1
            In_rot_yr+=1
            
        Abv_c_df_final=pd.DataFrame(merged_data['Year'], columns=['Year'])
    
        for i,j in zip (Abv_c_df.columns,rep):
            Abv_c_df_final[i]=j
            
        Abv_c_df_final=Abv_c_df_final.set_index('Year')
    
        c=0
        for j in Abv_c_df_final.columns:
            Abv_c_df_final[j]=Abv_c_df_final[j].shift(periods=c, fill_value=0)
            c+=1
        
        #############  Creating the final Dataframe ###########################
        Abv_c_df_final["Veg_C_Emissions"]=Abv_c_df_final.sum(axis=1)
    
        final_output_df= pd.DataFrame(Abv_c_df_final.index, columns=['Year'])
        final_output_df["LU_Area"]=np.array(merged_data["LU_Area_ha"])
        final_output_df["Agro_F_C_Seq"]=((np.array(Abv_c_df_final["Veg_C_Emissions"]))*3.67)/1000000
        final_output_df["Agro_S_C_Seq"]=((np.array(Soil_df["Soil_Emission"]))*3.67)/1000000
        final_output_df.loc[-1] = [df["Year"][0],df["LU_Area_ha"][0], 0,0]
        final_output_df.index = final_output_df.index + 1
        final_output_df = final_output_df.sort_index()
        final_output_df["Ini_F_VegC_Seq"]=((Ini_forest_C_Seq)*3.67)/1000000
        final_output_df["Ini_F_SoilC_Seq"]=((Ini_forest_Soil_C_Seq)*3.67)/1000000
        final_output_df["Total_Emissions"]=(final_output_df["Agro_F_C_Seq"]+final_output_df["Agro_S_C_Seq"]+final_output_df["Ini_F_VegC_Seq"]
                                            +final_output_df["Ini_F_SoilC_Seq"])
        
        label.config(text="Here is the Emissions/Sequestration Figure") 
        
        plt.bar(final_output_df['Year'],final_output_df["Total_Emissions"],color='green')
        plt.xlabel('Years')
        plt.ylabel('Emissions (MtCO2)')
        plt.ticklabel_format(style='plain')
        plt.show()
        
        label.config(text="Please save the output file")  
        file=asksaveasfilename(filetypes=[("csv file","*.csv")],defaultextension=".csv")
        final_output_df.to_csv(file, index=False)
        label.config(text="Processing over! The output file has been saved into the user-defined directory. Thanks!")
    

    else:
        Abv_c_df_final=pd.DataFrame(merged_data['Year'], columns=['Year'])

        Abv_c_df_final["Veg_C_Emissions"]=np.array(Abv_c_df.sum(axis=1))

        final_output_df= pd.DataFrame(merged_data['Year'], columns=['Year'])
        final_output_df["LU_Area"]=np.array(merged_data["LU_Area_ha"])
        final_output_df["Agro_F_C_Seq"]=((np.array(Abv_c_df_final["Veg_C_Emissions"]))*3.67)/1000000
        final_output_df["Agro_S_C_Seq"]=((np.array(Soil_df["Soil_Emission"]))*3.67)/1000000
        final_output_df.loc[-1] = [df["Year"][0],df["LU_Area_ha"][0], 0,0]
        final_output_df.index = final_output_df.index + 1
        final_output_df = final_output_df.sort_index()
        final_output_df["Ini_F_VegC_Seq"]=((Ini_forest_C_Seq)*3.67)/1000000
        final_output_df["Ini_F_SoilC_Seq"]=((Ini_forest_Soil_C_Seq)*3.67)/1000000
        
        final_output_df["Total_Emissions"]=(final_output_df["Agro_F_C_Seq"]+final_output_df["Agro_S_C_Seq"]+final_output_df["Ini_F_VegC_Seq"]
                                            +final_output_df["Ini_F_SoilC_Seq"])
        
        label.config(text="Here is the Emissions/Sequestration Figure") 
        
        plt.bar(final_output_df['Year'],final_output_df["Total_Emissions"],color='green')
        plt.xlabel('Years')
        plt.ylabel('Emissions (MtCO2)')
        plt.ticklabel_format(style='plain')
        plt.show()
        
        label.config(text="Please save the output file")  
        file=asksaveasfilename(filetypes=[("csv file","*.csv")],defaultextension=".csv")
        final_output_df.to_csv(file, index=False)
        label.config(text="Processing over! The output file has been saved into the user-defined directory. Thanks!")
    

####################  GUI Starts From Here.... ###############################################

window = Tk()
window.geometry("1000x600")
window.configure(bg = "#ffffff")
window.title("Agroforestry Emissions Calculator")
#window.wm_attributes("-transparentcolor",'blue')
canvas = Canvas(
    window,
    bg = "#ffffff",
    height = 600,
    width = 1000,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge")
canvas.place(x = 0, y = 0)

background_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\background.png")
background = canvas.create_image(
    500, 300,
    image=background_img)

entry0_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox0.png")
entry0_bg = canvas.create_image(
    650, 165,
    image = entry0_img)

entry0 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry0.place(
    x = 630, y = 144,
    width = 50.0,
    height = 42)

entry1_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox1.png")
entry1_bg = canvas.create_image(
    650, 220,
    image = entry1_img)

entry1 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry1.place(
    x = 630, y = 200,
    width = 50.0,
    height = 42)

entry2_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox2.png")
entry2_bg = canvas.create_image(
    650, 275,
    image = entry2_img)

entry2 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry2.place(
    x = 630, y = 252,
    width = 50.0,
    height = 45)

entry3_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox3.png")
entry3_bg = canvas.create_image(
    650, 330,
    image = entry3_img)

entry3 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry3.place(
    x = 630, y = 308,
    width = 50.0,
    height = 44)

entry4_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox4.png")
entry4_bg = canvas.create_image(
    650, 390,
    image = entry4_img)

entry4 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry4.place(
    x = 630, y = 365,
    width = 50.0,
    height = 49)

entry5_img = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img_textBox5.png")
entry5_bg = canvas.create_image(
    650, 450,
    image = entry5_img)

entry5 = Entry(
    bd = 0,
    bg = "#d9d9d9",
    highlightthickness = 0)

entry5.place(
    x = 630, y = 425,
    width = 50.0,
    height = 49)

img0 = PhotoImage(file = f"D:\\Agroforestry_App\\pyGUI\\Background_Files\\img0.png")
b0 = Button(
    image = img0,
    borderwidth = 0,
    highlightthickness = 0,
    command = btn_clicked,
    relief = "flat")

b0.place(
    x = 450, y = 500,
    width = 120,
    height = 37)
# Create a Label widget
label=Label(text="", font=('Calibri 15'))
#label.place(
#    x= 200 , y = 550,
#    width = 600,
#    height = 30
#)
label.pack()

window.resizable(False, False)
window.mainloop()
