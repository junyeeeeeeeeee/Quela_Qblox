import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import tkinter as tk
from tkinter import ttk
from Modularize.support import init_meas

QD_path = ''
dr = ''
ip = ''
mode = ''

# Initialize measurement window
def init_meas_window():
    
    root = tk.Tk()
    root.title('Init meas window')
    root.geometry('200x200')
    def init():
        global QD_path, dr, ip, mode
        QD_path = QD_path_en.get()
        dr = DR_box.get()
        ip = ip_box.get()
        mode = mode_box.get()
        if dr == '' or ip == '' or mode == '':
            raise KeyError('You have to insert the variables!')
        root.destroy()
    QD_path_la = tk.Label(root, text = 'QD path')
    QD_path_la.pack()
    QD_path_en = tk.Entry(root, width=30, justify='center')
    QD_path_en.pack()
    
    DR_la = tk.Label(root, text = 'DR location')   
    DR_la.pack()
    DR_box = ttk.Combobox(root,
                          width=15,
                          values=['DR1','DR2','DR3','DR4'])
    DR_box.pack()
    
    ip_la = tk.Label(root, text = 'IP location')   
    ip_la.pack()
    ip_box = ttk.Combobox(root,
                          width=15,
                          values=[170, 171])
    ip_box.pack()
    
    mode_la = tk.Label(root, text = 'mode')   
    mode_la.pack()
    mode_box = ttk.Combobox(root,
                          width=15,
                          values=['new','load'])
    mode_box.pack()

    btn = tk.Button(root, text='Enter', command=init)   # 建立按鈕，點擊按鈕時，執行 show 函式
    btn.pack()
    root.mainloop()
    
    Qmanager, cluster, meas_ctrl, ic = init_meas(QuantumDevice_path=QD_path,
                                                 dr_loc=dr,
                                                 cluster_ip=ip,
                                                 mode=mode)
    return Qmanager, cluster, meas_ctrl, ic, QD_path

