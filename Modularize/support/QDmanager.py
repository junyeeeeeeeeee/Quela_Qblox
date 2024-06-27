import os, datetime, pickle
from xarray import Dataset
from Modularize.support.FluxBiasDict import FluxBiasDict
from Modularize.support.Notebook import Notebook
from quantify_scheduler.device_under_test.quantum_device import QuantumDevice
from quantify_scheduler.device_under_test.transmon_element import BasicTransmonElement

def ret_q(dict_a):
    x = []
    for i in dict_a:
        if i[0] == 'q':
            x.append(i)
    return x

def ret_c(dict_a):
    x = []
    for i in dict_a:
        if i[0] == 'c':
            x.append(i)
    return x

class QDmanager():
    def __init__(self,QD_path:str=''):
        self.path = QD_path
        self.refIQ = {}
        self.Hcfg = {}
        self.Log = "" 
        self.Identity=""
        self.chip_name = ""
        self.chip_type = ""
        

    def register(self,cluster_ip_adress:str,which_dr:str,chip_name:str='',chip_type = ''):
        """
        Register this QDmanager according to the cluster ip and in which dr and the chip name.
        """
        specifier = cluster_ip_adress.split(".")[-1] # 192.168.1.specifier
        self.Identity = which_dr.upper()+"#"+specifier # Ex. DR2#171
        self.chip_name = chip_name
        self.chip_type = chip_type

    def memo_refIQ(self,ref_dict:dict):
        """
        Memorize the reference IQ according to the given ref_dict, which the key named in "q0"..., and the value is composed in list by IQ values.\n
        Ex. ref_dict={"q0":[0.03,-0.004],...} 
        """
        for q in ref_dict:
            self.refIQ[q] = ref_dict[q]
    
    def refresh_log(self,message:str):
        """
        Leave the message for this file.
        """
        self.Log = message

    def QD_loader(self, new_Hcfg:bool=False):
        """
        Load the QuantumDevice, Bias config, hardware config and Flux control callable dict from a given json file path contain the serialized QD.
        """
        with open(self.path, 'rb') as inp:
            gift = pickle.load(inp) # refer to `merged_file` in QD_keeper()
        # string and int
        self.chip_name:str = gift["chip_info"]["name"]
        self.chip_type:str = gift["chip_info"]["type"]
        self.Identity:str = gift["ID"]
        self.Log:str = gift["Log"]
        self.q_num:int = len(list(filter(ret_q,gift["Flux"])))
        self.c_num:int = len(list(filter(ret_c,gift["Flux"])))
        # class    
        self.Fluxmanager :FluxBiasDict = FluxBiasDict(qb_number=self.q_num,cp_number=self.c_num)
        self.Fluxmanager.activate_from_dict(gift["Flux"])
        self.Notewriter: Notebook = Notebook(q_number=self.q_num)
        self.Notewriter.activate_from_dict(gift["Note"])
        self.quantum_device :QuantumDevice = gift["QD"]
        # dict
        if new_Hcfg:
            from Modularize.support.Experiment_setup import hcfg_map
            self.Hcfg = hcfg_map[self.Identity.split("#")[0].lower()]
        else:
            self.Hcfg = gift["Hcfg"]
        self.refIQ:dict = gift["refIQ"]
        
        self.quantum_device.hardware_config(self.Hcfg)
        print("Old friends loaded!")
    
    def QD_keeper(self, special_path:str=''):
        """
        Save the merged dictionary to a json file with the given path. \n
        Ex. merged_file = {"QD":self.quantum_device,"Flux":self.Fluxmanager.get_bias_dict(),"Hcfg":Hcfg,"refIQ":self.refIQ,"Log":self.Log}
        """
        if self.path == '' or self.path.split("/")[-2].split("_")[-1] != datetime.datetime.now().day:
            db = Data_manager()
            db.build_folder_today()
            self.path = os.path.join(db.raw_folder,f"{self.Identity}_SumInfo.pkl")
        Hcfg = self.quantum_device.generate_hardware_config()
        # TODO: Here is onlu for the hightlighs :)
        merged_file = {"ID":self.Identity,"chip_info":{"name":self.chip_name,"type":self.chip_type},"QD":self.quantum_device,"Flux":self.Fluxmanager.get_bias_dict(),"Hcfg":Hcfg,"refIQ":self.refIQ,"Note":self.Notewriter.get_notebook(),"Log":self.Log}
        
        with open(self.path if special_path == '' else special_path, 'wb') as file:
            pickle.dump(merged_file, file)
            print(f'Summarized info had successfully saved to the given path!')

    

    def build_new_QD(self,qubit_number:int,coupler_number:int,Hcfg:dict,cluster_ip:str,dr_loc:str,chip_name:str='',chip_type:str=''):

        """
        Build up a new Quantum Device, here are something must be given about it:\n
        (1) qubit_number: how many qubits is in the chip.\n
        (2) Hcfg: the hardware configuration between chip and cluster.\n
        (3) cluster_ip: which cluster is connected. Ex, cluster_ip='192.168.1.171'\n
        (4) dr_loc: which dr is this chip installed. Ex, dr_loc='dr4'
        """
        print("Building up a new quantum device system....")
        self.q_num = qubit_number
        self.cp_num = coupler_number
        self.Hcfg = Hcfg
        self.chip_name = chip_name
        self.chip_type = chip_type
        self.register(cluster_ip_adress=cluster_ip,which_dr=dr_loc,chip_name=chip_name,chip_type=chip_type)

        self.Fluxmanager :FluxBiasDict = FluxBiasDict(self.q_num,self.cp_num)
        self.Notewriter: Notebook = Notebook(self.q_num)
        """ #for firmware v0.6.2
        self.quantum_device = QuantumDevice("academia_sinica_device")
        self.quantum_device.hardware_config(self.Hcfg)
        
        # store references
        self.quantum_device._device_elements = list()

        for i in range(qubit_number):
            qubit = BasicTransmonElement(f"q{i}")
            qubit.measure.acq_channel(i)
            self.quantum_device.add_element(qubit)
            self.quantum_device._device_elements.append(qubit)
        """
        # for firmware v0.7.0
        from qcodes.instrument import find_or_create_instrument
        self.quantum_device = find_or_create_instrument(QuantumDevice, recreate=True, name=f"QPU{dr_loc.lower()}")
        self.quantum_device.hardware_config(self.Hcfg)
        for qb_idx in range(self.q_num):
            qubit = find_or_create_instrument(BasicTransmonElement, recreate=True, name=f"q{qb_idx}")
            qubit.measure.acq_channel(qb_idx)
            self.quantum_device.add_element(qubit)
        
    def keep_meas_option(self,target_q:str,z_bias:float,modi_idx:int):
        """ keep the following info into Notebook\n
        1) XY freq.\n
        2) RO freq.\n
        3) RO amp.\n
        4) pi-pulse amp.\n
        5) 2tone_pi amp.\n
        6) pi-pulse duration.\n
        7) ref-IQ point.\n
        8) bias of this point.\n
        9) ro attenuation.
        """
        print(z_bias)
        if modi_idx != "-1":
            if len(self.Notewriter.get_all_meas_options(target_q)) <= modi_idx:
                self.Notewriter.create_meas_options(target_q)
        qubit = self.quantum_device.get_element(target_q)
        ROF = qubit.clock_freqs.readout()
        XYF = qubit.clock_freqs.f01()
        pi_amp = qubit.rxy.amp180()
        conti_pi_amp = self.Notewriter.get_2tone_piampFor(target_q)
        pi_dura = qubit.rxy.duration()
        ref_iq = self.refIQ[target_q]
        ro_attenuation = self.Notewriter.get_DigiAtteFor(target_q,'ro')
        ro_amp = qubit.measure.pulse_amp()
        option_dict = {"f01":XYF,"rof":ROF,"rop":ro_amp,"pi_amp":pi_amp,"2tone_pi_amp":conti_pi_amp,"pi_dura":pi_dura,"refIQ":ref_iq,"bias":z_bias,"ro_atte":ro_attenuation}

        self.Notewriter.write_meas_options({target_q:option_dict},modi_idx)
        print(f"Optional meas point had been recorded! @ Z~{round(z_bias,3)}")


    def write_with_meas_option(self,target_q:str,idx_chosen:str):
        """ call the following info into QuantumDevice, Fluxmanager, Notewriter, QDmanager\n
        1) XY freq.\n
        2) RO freq.\n
        3) RO amp.\n
        4) pi-pulse amp.\n
        5) 2tone_pi amp.\n
        6) pi-pulse duration.\n
        7) ref-IQ point.\n
        8) bias of this point.\n
        9) ro attenuation.
        """
        option_selected = self.Notewriter.get_all_meas_options(target_q)[int(idx_chosen)]
        qubit = self.quantum_device.get_element(target_q)
        qubit.clock_freqs.readout(option_selected["rof"])
        qubit.clock_freqs.f01(option_selected["f01"])
        qubit.rxy.amp180(option_selected["pi_amp"])
        self.Notewriter.save_2tone_piamp_for(target_q,float(option_selected["2tone_pi_amp"]))
        qubit.rxy.duration(option_selected["pi_dura"])
        self.refIQ[target_q] = option_selected["refIQ"]
        self.Notewriter.save_DigiAtte_For(int(option_selected["ro_atte"]),target_q,'ro')
        qubit.measure.pulse_amp(option_selected["rop"])
        if idx_chosen != '0':
            self.Fluxmanager.save_tuneawayBias_for('manual',target_q,option_selected["bias"])
            self.Fluxmanager.press_offsweetspot_button(target_q,True) # here is the only way to press this button
        else:
            self.Fluxmanager.save_sweetspotBias_for(target_q,option_selected["bias"])
            self.Fluxmanager.press_offsweetspot_button(target_q,False)

    ### Convenient short cuts
# Object to manage data and pictures store.

class Data_manager:
    
    def __init__(self):
        from Modularize.support.Path_Book import meas_raw_dir
        from Modularize.support.Path_Book import qdevice_backup_dir
        self.QD_back_dir = qdevice_backup_dir
        self.raw_data_dir = meas_raw_dir

    # generate time label for netCDF file name
    def get_time_now(self)->str:
        """
        Since we save the Xarray into netCDF, we use the current time to encode the file name.\n
        Ex: 19:23:34 return H19M23S34 
        """
        current_time = datetime.datetime.now()
        return f"H{current_time.hour}M{current_time.minute}S{current_time.second}"
    
    def get_date_today(self)->str:
        current_time = datetime.datetime.now()
        return f"{current_time.year}_{current_time.month}_{current_time.day}"

    # build the folder for the data today
    def build_folder_today(self,parent_path:str=''):
        """
        Build up and return the folder named by the current date in the parent path.\n
        Ex. parent_path='D:/Examples/'
        """ 
        if parent_path == '':
            parent_path = self.QD_back_dir

        folder = self.get_date_today()
        new_folder = os.path.join(parent_path, folder) 
        if not os.path.isdir(new_folder):
            os.mkdir(new_folder) 
            print(f"Folder {folder} had been created!")

        pic_folder = os.path.join(new_folder, "pic")
        if not os.path.isdir(pic_folder):
            os.mkdir(pic_folder) 
        
        self.raw_folder = new_folder
        self.pic_folder = pic_folder

    
    def save_raw_data(self,QD_agent:QDmanager,ds:Dataset,qb:str='q0',label:str=0,exp_type:str='CS', specific_dataFolder:str='', get_data_loc:bool=False):
        exp_timeLabel = self.get_time_now()
        self.build_folder_today(self.raw_data_dir)
        parent_dir = self.raw_folder if specific_dataFolder =='' else specific_dataFolder
        dr_loc = QD_agent.Identity.split("#")[0]
        if exp_type.lower() == 'cs':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_CavitySpectro_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'pd':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_PowerCavity_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'fd':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_FluxCavity_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'ss':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_SingleShot({label})_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == '2tone':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_2tone_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'f2tone':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_Flux2tone_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'powerrabi':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_powerRabi_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'timerabi':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_timeRabi_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'ramsey':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_ramsey_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 't1':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_T1({label})_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 't2':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_T2({label})_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'rofcali':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_RofCali({label})_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        elif exp_type.lower() == 'zt1':
            path = os.path.join(parent_dir,f"{dr_loc}{qb}_zT1({label})_{exp_timeLabel}.nc")
            ds.to_netcdf(path)
        else:
            path = ''
            raise KeyError("Wrong experience type!")
        
        if get_data_loc:
            return path
    
    def save_histo_pic(self,QD_agent:QDmanager,hist_dict:dict,qb:str='q0',mode:str="t1", show_fig:bool=False, save_fig:bool=True,pic_folder:str=''):
        from Modularize.support.Pulse_schedule_library import hist_plot
        exp_timeLabel = self.get_time_now()
        if pic_folder == '':
            self.build_folder_today(self.raw_data_dir)
            pic_dir = self.pic_folder
        else:
            pic_dir = pic_folder
        dr_loc = QD_agent.Identity.split("#")[0]
        if mode.lower() =="t1" :
            if save_fig:
                fig_path = os.path.join(pic_dir,f"{dr_loc}{qb}_T1histo_{exp_timeLabel}.png")
            else:
                fig_path = ''
            hist_plot(qb,hist_dict ,title=f"T1",save_path=fig_path, show=show_fig)
        elif mode.lower() =="t2" :
            if save_fig:
                fig_path = os.path.join(pic_dir,f"{dr_loc}{qb}_T2histo_{exp_timeLabel}.png")
            else:
                fig_path = ''
            hist_plot(qb,hist_dict ,title=f"T2",save_path=fig_path, show=show_fig)
        elif mode.lower() in ["ss", "os"] :
            if save_fig:
                fig_path = os.path.join(pic_dir,f"{dr_loc}{qb}_effThisto_{exp_timeLabel}.png")
            else:
                fig_path = ''
            hist_plot(qb,hist_dict ,title=f"eff_T",save_path=fig_path, show=show_fig)
        else:
            raise KeyError("mode should be 'T1' or 'T2'!")
        
    def save_dict2json(self,QD_agent:QDmanager,data_dict:dict,qb:str='q0',get_json:bool=False):
        """
        Save a dict into json file. Currently ONLY support z-gate 2tone fitting data.
        """
        import json
        exp_timeLabel = self.get_time_now()
        self.build_folder_today(self.raw_data_dir)
        dr_loc = QD_agent.Identity.split("#")[0]
        path = os.path.join(self.raw_folder,f"{dr_loc}{qb}_FluxFqFIT_{exp_timeLabel}.json")
        with open(path, "w") as json_file:
            json.dump(data_dict, json_file)
        print("Flux vs fq to-fit data had been saved!")
        if get_json:
            return path
    
    def get_today_picFolder(self)->str:
        """
        Get the picture folder today. Return its path.
        """
        self.build_folder_today(self.raw_data_dir)
        return self.pic_folder


