def URF():
    depthmyr = pd.read_csv("E:/")
    depthmyr = depthmyr.values.tolist()
    outputpath = "E:/"
    reach = pd.read_csv("E:/")
    reach = reach.values.tolist()
    for i in range(3, len(depthmyr) - 3, 1):

        
        if depthmyr[i-3][32] != -100 and depthmyr[i+3][32] != -100:
            OcLUp = STRtoLIST(depthmyr[i-3][38])
            OcLDn = STRtoLIST(depthmyr[i+3][38])
            PSUp = depthmyr[i-3][39]
            PSDn = depthmyr[i+3][39]
            
            #OcLUp = [x for x in OcLUp if x != 0]
            #OcLDn = [x for x in OcLDn if x != 0]
            
            frequency_counter_Up = Counter(OcLUp)
            total_elements_Up = len(OcLUp)
            frequency_percentage_Up = [[key, float(count)/total_elements_Up] for key, count in frequency_counter_Up.items()]
            
            frequency_counter_Dn = Counter(OcLDn)
            total_elements_Dn = len(OcLDn)
            frequency_percentage_Dn = [[key, float(count)/total_elements_Dn] for key, count in frequency_counter_Dn.items()]
            
            Slope = []
            Freqc = []
            #print(OcLUp)
            #print(frequency_percentage_Up)
            for j in range(len(frequency_percentage_Up)):
                
                if frequency_percentage_Up[j][0] == 0:
                    
                    LUp = depthmyr[i-3][34]
                    FreUp = frequency_percentage_Up[j][1]
                    
                else:
                
                    lnWUp = math.log(PSUp * frequency_percentage_Up[j][0])
                    FreUp = frequency_percentage_Up[j][1]
                 
                    lnDUp =  (lnWUp * depthmyr[i-3][32]) + depthmyr[i-3][33]
                
                    LUp = math.exp(lnDUp) + depthmyr[i-3][34]
                
                for u in range(len(frequency_percentage_Dn)):
                    
                    if PSDn * frequency_percentage_Dn[u][0] == 0:
                        
                        LDn = depthmyr[i+3][34]
                        FreDn = frequency_percentage_Dn[u][1]
                        
                    else:
                    
                        lnWDn = math.log(PSDn * frequency_percentage_Dn[u][0])
                        FreDn = frequency_percentage_Dn[u][1]
                    
                        lnDDn =  (lnWDn * depthmyr[i+3][32]) + depthmyr[i+3][33]

                        LDn = math.exp(lnDDn) + depthmyr[i+3][34]
                    #print(LUp, LDn)
                    Dis = (depthmyr[i-3][19] / 2.0) + depthmyr[i-2][19] + depthmyr[i-1][19] + depthmyr[i][19] + depthmyr[i+1][19] + depthmyr[i+2][19] + (depthmyr[i+3][19] / 2.0)
                    
                    S = (LUp - LDn) / Dis
                    if S < 0.0:
                        S = -1.0
                    Slope.append(S)
                    Freqc.append(FreUp * FreDn)
                    
                
                
            #print(Dis)
            #print(Slope)
            #print(Freqc)
            #print(sum(Freqc))  
            combined_lists = [[x, y] for x, y in zip(Slope, Freqc)]
            newslopelist = [sublist for sublist in combined_lists if sublist[0] != -1]
            
            if len(newslopelist) != 0:
                newpro = sum(sublist[1] for sublist in newslopelist)
                newcombined_lists = [[sublist[0], float(sublist[1] / newpro)] for sublist in newslopelist]
                #print(newcombined_lists)
                for k in range(len(reach)):
                    if depthmyr[i][20] == int(round(reach[k][3])):
                        true_value = float(reach[k][16]/1000.0)
                        break;

                squared_errors_sum = 0
                for measurement in newcombined_lists:
                    value, probability = measurement
                    error = value - true_value
                    squared_errors_sum += (error ** 2) * probability
                
                rmse = math.sqrt(squared_errors_sum)
            
                E_list = [a * b for a, b in newcombined_lists]
                E = sum(E_list)
                depthmyr[i][40] = E
                depthmyr[i][41] = rmse
                print(E)
            
            else:
                depthmyr[i][40] = -100
                depthmyr[i][41] = -100
            
            
        else:
            depthmyr[i][40] = -100
            depthmyr[i][41] = -100
            continue;   

    filtered_matrix = [row for row in depthmyr if all(row[i] != -100 for i in [30, 32, 40])]
    test = pd.DataFrame(data = filtered_matrix)
    test.to_csv(outputpath + "", index = False, header = True)

def Widthoc1(OcList, stdvalue):
    Fi = []
    R = []
    
    Wmax = float(OcList[0])
    for k in range(0, len(OcList), 1):
        if OcList[k] != 0:
            Wi = float(OcList[k])
            Fi.append(float(Wi/Wmax)*100.0)
            #R.append(-math.log((k*5.0)/100.0))
    
    Fi = np.array(Fi)
    stdvalue = np.array(stdvalue)
    
    rmse = np.sqrt(np.mean((stdvalue - Fi)**2))   

    return rmse

def FSE():
    depth = pd.read_csv("E:/")
    outputpath = "E:/"
    depth = depth.values.tolist()
    for i in range(3, len(depth)-3, 1):
        
        OcL1 = STRtoLIST(depth[i-3][35])
        OcL2 = STRtoLIST(depth[i-2][35])
        OcL3 = STRtoLIST(depth[i-1][35])
        OcL4 = STRtoLIST(depth[i][35])
        OcL5 = STRtoLIST(depth[i+1][35])
        OcL6 = STRtoLIST(depth[i+2][35])
        OcL7 = STRtoLIST(depth[i+3][35])
        
        combined_lists = zip(OcL1, OcL2, OcL3, OcL4, OcL5, OcL6, OcL7)
        
        OcL = [sum(x) for x in combined_lists]
        
        #OcL = STRtoLIST(depth[i][10])#16
        if OcL[0] == 0:
            depth[i][36] = -100000000
            depth[i][37] = -100000000
            continue;
        n0 = 0
        for j in range(0, len(OcL), 1):
            if OcL[j] == 0.0:#if OcL[j] / OcL[0] <= 0.1:
                n0 = n0 + 1
        if n0 == 0:
            stdvalue = [100.0, 95.0, 90.0, 85.0, 80.0, 75.0, 70.0, 65.0, 60.0, 55.0, 50.0, 45.0, 40.0, 35.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0]
        elif n0 == 1:
            stdvalue = [100.0, 94.737, 89.474, 84.211, 78.947, 73.684, 68.421, 63.158, 57.895, 52.632, 47.368, 42.105, 36.842, 31.579, 26.316, 21.053, 15.789, 10.526, 5.263]
        elif n0 == 2:
            stdvalue = [100.0, 94.444, 88.889, 83.333, 77.778, 72.222, 66.667, 61.111, 55.556, 50.0, 44.444, 38.889, 33.333, 27.778, 22.222, 16.667, 11.111, 5.556]
        elif n0 == 5:
            stdvalue = [100.0, 93.333, 86.667, 80.0, 73.333, 66.667, 60.0, 53.333, 46.667, 40.0, 33.333, 26.667, 20.0, 13.333, 6.667]
        elif n0 == 4:
            stdvalue = [100.0, 93.75, 87.5, 81.25, 75.0, 68.75, 62.5, 56.25, 50.0, 43.75, 37.5, 31.25, 25.0, 18.75, 12.5, 6.25]
        elif n0 == 3:
            stdvalue = [100.0, 94.118, 88.235, 82.353, 76.471, 70.588, 64.706, 58.824, 52.941, 47.059, 41.176, 35.294, 29.412, 23.529, 17.647, 11.765, 5.882]
        elif n0 == 8:
            stdvalue = [100.0, 91.667, 83.333, 75.0, 66.667, 58.333, 50.0, 41.667, 33.333, 25.0, 16.667, 8.333]
        elif n0 == 9:
            stdvalue = [100.0, 90.909, 81.818, 72.727, 63.636, 54.545, 45.455, 36.364, 27.273, 18.182, 9.091]
        elif n0 == 10:
            stdvalue = [100.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0]
        elif n0 == 11:
            stdvalue = [100.0, 88.889, 77.778, 66.667, 55.556, 44.444, 33.333, 22.222, 11.111]
        elif n0 == 7:
            stdvalue = [100.0, 92.308, 84.615, 76.923, 69.231, 61.538, 53.846, 46.154, 38.462, 30.769, 23.077, 15.385, 7.692]
        elif n0 == 6:
            stdvalue = [100.0, 92.857, 85.714, 78.571, 71.429, 64.286, 57.143, 50.0, 42.857, 35.714, 28.571, 21.429, 14.286, 7.143]
        elif n0 == 12:
            stdvalue = [100.0, 87.5, 75.0, 62.5, 50.0, 37.5, 25.0, 12.5]
        elif n0 == 13:
            stdvalue = [100.0, 85.714, 71.429, 57.143, 42.857, 28.571, 14.286]
        elif n0 == 14:
            stdvalue = [100.0, 83.333, 66.667, 50.0, 33.333, 16.667]
        elif n0 == 15:
            stdvalue = [100.0, 80.0, 60.0, 40.0, 20.0]
        elif n0 == 16:
            stdvalue = [100.0, 75.0, 50.0, 25.0]
        elif n0 == 17:
            stdvalue = [100.0, 66.667, 33.333]
        else:
            depth[i][36] = -100000000
            depth[i][37] = -100000000
            continue;   
        WF1 = Widthoc1(OcL, stdvalue)
        depth[i][36] = WF1
        depth[i][37] = WF1
        print('done', i)
       
        
    #depth.sort(key=lambda x:x[7], reverse=False)
    test = pd.DataFrame(data = depth)
    test.to_csv(outputpath + "", index = False, header = True)

def SCT():
    SCT = pd.read_csv("D:/")
    SCT = SCT.values.tolist()
    outputpath = "D:/"
    
    
    for i in range(3, len(SCT)-3, 1):#3, len(SCT)-3, 1
        K1 = []
        B1 = []
        for j in range(-3, 4, 1):
            #print(i+j)
            Cross_sect = STRtoLIST(SCT[i+j][4])
            wse = min(Cross_sect)
            PS = 700.0 / (len(Cross_sect) - 1)
            if (wse * 10) % 5 == 0:
                if 1==0:#Cross_sect[0] == wse or Cross_sect[len(Cross_sect) - 1] == wse:
                    #print('1')
                    continue;
                else:
                    W = []
                    L = []
                    io_list = []
                    ioo_list = []
                    order1_list = []
                    order2_list = []
                    if Cross_sect[0] == wse:
                        order1_list.append(0)
                        
                    for k1 in range(1, len(Cross_sect), 1):
                        if Cross_sect[k1] == wse and Cross_sect[k1-1] != wse:
                            order1_list.append(k1)
                            #break;
                        
                    for k2 in range(0, len(Cross_sect) - 1, 1):
                        if Cross_sect[k2] == wse and Cross_sect[k2+1] != wse:
                            order2_list.append(k2)
                            #break;
                            
                    if Cross_sect[len(Cross_sect) - 1] == wse:
                        order2_list.append(len(Cross_sect) - 1)
                    
                    order_dis = [x - y for x, y in zip(order2_list, order1_list)]
                    
                    order_index = order_dis.index(max(order_dis))
                    
                    order1 = order1_list[order_index]
                    order2 = order2_list[order_index]
                    
                    if order2 - order1 == 0:
                        SCT[i+j][32] = -100
                        SCT[i+j][33] = -100
                        SCT[i+j][34] = -100
                        continue;
                    else:
                        
                        
                        #print(order1, order2)
                        #W.append((order2 - order1) * 25.0)
                        #L.append(wse)
                        
                        for io1 in range(order1 - 1, 0, -1):
                            
                            tr = Cross_sect[io1]
                            st = -1
                            ioo = -1
                            
                            for io1_1 in range(io1 + 1, len(Cross_sect) - 1, 1):
                                
                                if Cross_sect[io1_1] >= tr:
                                    
                                    st = Cross_sect[io1_1]
                                    ioo = io1_1
                                    break;
                            
                            if ioo == -1 or ioo < order2:
                                ioo_list.append(-1)
                                io_list.append(-1)
                                W.append(0)
                                L.append(0)
                                
                            else:
                                #PS = 700.0 / (len(Cross_sect) - 1)
                                ioo_list.append(ioo)
                                io_list.append(io1)
                                W.append((tr - Cross_sect[ioo - 1]) / (st - Cross_sect[ioo - 1]) * PS + ((ioo - io1 - 1) * PS))
                                L.append(tr)
                                
                        for io2 in range(order2 + 1, len(Cross_sect), 1):
                            
                            tr1 = Cross_sect[io2]
                            st1 = -1
                            ioo1 = -1
                            
                            for io2_1 in range(io2 - 1, 0, -1):
                                
                                if Cross_sect[io2_1] >= tr1:
                                    
                                    st1 = Cross_sect[io2_1]
                                    ioo1 = io2_1
                                    break;
                            
                            if ioo1 == -1 or ioo1 > order1:
                                
                                ioo_list.append(-1)
                                io_list.append(-1)
                                W.append(0)
                                L.append(0)
                                
                            else:
                                
                                #PS = 700.0 / (len(Cross_sect) - 1)
                                ioo_list.append(ioo1)
                                io_list.append(io2)
                                W.append((tr1 - Cross_sect[ioo1 + 1]) / (st1 - Cross_sect[ioo1 + 1]) * PS + ((io2 - ioo1 - 1) * PS))
                                L.append(tr1)
                            
                                #print(botriver)
            #W = W.sort()
            #L = L.sort()
            #print(W)
            #print(L)
            #print(sorted(W))
            #print(sorted(L))
            else:
                W = [0]
                L = [0]
                io_list = [0]
                ioo_list = [0]
            
            if all(x == 0 for x in W):
                
                SCT[i+j][32] = -100
                SCT[i+j][33] = -100
                SCT[i+j][34] = -100
                continue;

            two_d_list = [[w, x, y, z] for w, x, y, z in zip(io_list, W, L, ioo_list)]
            sorted_two_d_list = sorted(two_d_list, key=lambda x: x[1])
            
            BankFullW = SCT[i+j][28]
            BankFullD = SCT[i+j][29]
            #PS = 700.0 / (len(Cross_sect) - 1)
            
            if BankFullW <= (order2 - order1) * PS:
                
                Bathem = wse - BankFullD
                
            elif BankFullW >= sorted_two_d_list[-1][1]:
                
                Bathem = sorted_two_d_list[-1][2] - BankFullD
                
            else:
                
                for u in range(len(sorted_two_d_list)):
                    #print(BankFullW, sorted_two_d_list[u][1])
                    if BankFullW < sorted_two_d_list[u][1]:
                        #bfio = u
                        #print(u)
                        break;
                bfio = u    
                #print(len(sorted_two_d_list))
                #PS = 700.0 / (len(Cross_sect) - 1)
                line1_index1 = min(sorted_two_d_list[bfio][0], sorted_two_d_list[bfio][3])
                line2_index2 = max(sorted_two_d_list[bfio][0], sorted_two_d_list[bfio][3])
                line1_index2 = line1_index1 + 1
                line2_index1 = line2_index2 - 1
                
                line1_k = (Cross_sect[line1_index2] - Cross_sect[line1_index1])/PS
                line1_b = Cross_sect[line1_index1] - (line1_k * (line1_index1 * PS))
                
                line2_k = (Cross_sect[line2_index2] - Cross_sect[line2_index1])/PS
                line2_b = Cross_sect[line2_index1] - (line2_k * (line2_index1 * PS))
                
                BankFullL = (BankFullW*line1_k*line2_k - (line1_b*line2_k - (line2_b*line1_k))) / (line1_k - line2_k)
                    
                if BankFullL < sorted_two_d_list[bfio - 1][2]:
                    BankFullL = sorted_two_d_list[bfio - 1][2]
                   
                Bathem = BankFullL - BankFullD
                
            #print(Bathem)
            
            W.append((order2 - order1) * PS)
            L.append(wse)
            
            filtered_W = [value for value in sorted(W) if value != 0]
            filtered_L = [value for value in sorted(L) if value != 0]
            
            filtered_D = [x - Bathem for x in filtered_L]
            #print(filtered_W, filtered_L)
            filtered_W1 = []
            filtered_D1 = []
            
            for e in range(len(filtered_D)):
                
                if filtered_D[e] > 0.0:
                    filtered_D1.append(filtered_D[e])
                    filtered_W1.append(filtered_W[e])
            
            # Generate some example data
            
            if len(filtered_D1) < 3:
                SCT[i+j][32] = -100
                SCT[i+j][33] = -100
                SCT[i+j][34] = -100
                continue;
                
            W_data1 = [math.log(x) for x in filtered_W1]
            D_data1 = [math.log(x) for x in filtered_D1]
            
            #print(filtered_W1)
            #print(filtered_D1)
            
            W_data = np.array(W_data1)
            D_data = np.array(D_data1)
            
            k, b = np.polyfit(W_data, D_data, 1)
            SCT[i+j][32] = k
            SCT[i+j][33] = b
            SCT[i+j][34] = Bathem
            B1.append(b)
            K1.append(k)
        if len(B1) > 4:
            if SCT[i][32] == -100:
                SCT[i][30] = -100
                SCT[i][31] = -100
            else:
                E_std = SCT[i][32]
                C_std = math.exp(SCT[i][33])
                AVE = []
            
                for v in range(-3, 4, 1):
                    if v == 0 or SCT[i+v][32] == -100:
                        continue;
                
                    x = np.linspace(0.0001, 700, int(700/PS))
                    y_std = C_std * (x**E_std)
                    y_comp = math.exp(SCT[i+v][33]) * (x**SCT[i+v][32])
                    
                    
                    y_std = np.array(y_std)
                    y_comp = np.array(y_comp)
                    
                    rmse = np.sqrt(np.mean((y_comp - y_std)**2))
                    
                    rrmse = rmse / np.mean(y_std)
                    
                    #print(SCT[i+v][32])
                    #relative_error = np.abs(y_std - y_comp) / np.maximum(y_std, y_comp)
                    #average_relative_error = np.mean(relative_error)
                    AVE.append(rrmse)
                
                #print(AVE)
                SCT[i][30] = np.mean(AVE)
                SCT[i][31] = np.mean(AVE)
        else:
            SCT[i][30] = -100
            SCT[i][31] = -100
        #print(W_data)
        #print(L_data)
        #print(A_fit, k_fit, B_fit)
        print(i, SCT[i][0])
    test = pd.DataFrame(data = SCT)
    test.to_csv(outputpath + "", index = False, header = True)

def MSI():
    depth = pd.read_csv("E:/")
    outputpath = "E:/"
    depth = depth.values.tolist()
    Rw = []
    S = []
    LT = []
    LA = []
    LA1 = []
    for m in range(0, len(depth)):#len(depth)
        Rw.append(abs(depth[m][36]))#2.520681
        S.append(depth[m][41])
        LT.append(depth[m][30])
        LA.append(depth[m][31])
        #LA1.append(depth[m][2])
    
    minRW = (np.min(Rw) - np.mean(Rw)) / np.std(Rw)
    minS = (np.min(S) - np.mean(S)) / np.std(S)
    minLT = (np.min(LT) - np.mean(LT)) / np.std(LT)
    minLA = (np.min(LA) - np.mean(LA)) / np.std(LA)
    
    
    
    for i in range(0, len(depth)):
        d1 = (abs(depth[i][36]) - np.mean(Rw)) / np.std(Rw) - minRW
        d2 = (depth[i][41] - np.mean(S)) / np.std(S) - minS
        d3 = (depth[i][30] - np.mean(LT)) / np.std(LT) - minLT
        d4 = (depth[i][31] - np.mean(LA)) / np.std(LA) - minLA
        #d5 = (depth[i][2] - np.mean(LA1)) / np.std(LA1)
        ##if d1 < 0.0 and d2 < 0.0 and d3 < 0.0 and d4 < 0.0:
        #print(d3)
        depth[i][42] = math.sqrt((0.453 * d1) ** 2 + ((0.421 * d3) ** 2) + ((0.126 * d2) ** 2))
        print(depth[i][0], d1,d3,d2)
    depth.sort(key=lambda x:x[42], reverse = False)
    test = pd.DataFrame(data = depth)
    test.to_csv(outputpath + "", index = False, header = True)
