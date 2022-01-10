import preprocessData
import calculate
from calculate import Math
import numpy as np
import database
import sys
import os
import csv

class electrostatic:
    def __init__(self, frag2D, fragStruct3D, basis):
        path = './componentDatabase/ele_' + basis + '.csv'
        with open(path, newline='') as f:
            rows = csv.reader(f, delimiter=',')
            eleParameter = np.asarray(list(rows))
        self.path = path
        self.frag = frag2D
        self.fragStruct = fragStruct3D
        self.parameterEle = eleParameter
    def signEle(self):
        sign2d = []
        AB = ['A', 'B']
        CD = ['C', 'D']
        flag = True
        cdflag = False
        for frag1d in self.frag:
            for i, f in enumerate(frag1d):
                if f not in AB:
                    if f in CD:
                        startind = i
                        cdflag = True
                        break
                    startind = i
                    flag = False
                    break
            if flag:
                if cdflag:
                    PN = 1
                else:
                    startind = 0
                    PN = 1
            else:
                if len(sign2d) == 0:
                    PN = 1
                else:
                    if len(frag1d) == 1 and frag1d[0] == 'A':
                        startind = 0
                    PN = -1
            l, r = 0, len(frag1d)-1
#             print(l, r)
            lflag, rflag = startind, startind
            sign = [PN]
            while True:
                PN *= -1
                rflag += 1
                lflag -= 1
#                 print(rflag, lflag)
                if rflag > r and lflag < l:
                    break
                if rflag <= r:
                    sign.append(PN)
                if lflag >= l:
                    sign.insert(0, PN)
            sign2d.append(sign)
        print(sign2d)
        print(self.frag)
        self.sign2d = sign2d
    
    def energyEle(self, Rphi, Rrange, beta, train_mode, referenceEle):
        component, count_component, dis_component = [], [], []
        for i, i_frag in enumerate(self.fragStruct[0]):
            for j, j_frag in enumerate(self.fragStruct[1]):
                if self.sign2d[0][i] * self.sign2d[1][j] < 0:
                    dis = Math.distance3D(i_frag, j_frag)
                    err = Math.absoluteErr(dis, Rphi)
                    if err < Rrange:
                        if self.frag[0][i] < self.frag[1][j]:
                            com = self.frag[0][i] + self.frag[1][j]
                        else:
                            com = self.frag[1][j] + self.frag[0][i]
                        if com not in component:
                            component.append(com)
                            count_component.append(1)
                            dis_component.append([dis])
                        else:
                            com_ind = component.index(com)
                            count_component[com_ind] += 1
                            dis_component[com_ind].append(dis)
        print('Ele_component: ', component, '\n', 
              'Ele_component_quantity: ', count_component, '\n',
              'Ele_component_distance: ', dis_component, '\n')
        
        if train_mode:
            self.component = component
            self.count_component = count_component
            self.dis_component = dis_component
            print(self.modeling(Rphi, Rrange, beta, referenceEle))
            return self.parameterEle
        
        
        total = 0
        for i, comp in enumerate(component):
            ind = self.parameterEle[:,0].tolist().index(comp)
            energy = float(self.parameterEle[:,1][ind])
            distance = float(self.parameterEle[:,2][ind])
            beta1 = float(self.parameterEle[:,3][ind])
            for dis in dis_component[i]:
                e = energy*(2-(dis/distance)**beta1)
                total += e
        return total
    
    def modeling(self, Rphi, Rrange, beta, referenceEle):
        try:    
            complist = self.parameterEle[:,0].tolist() 
        except:
            complist = []
        unknow, dis_unknow = [], []
        for i, comp in enumerate(self.component):
            if comp not in complist:
                unknow.append(comp)
                dis_unknow.append(self.dis_component[i])
            else:
                ind = self.parameterEle[:,0].tolist().index(comp)
                energy = float(self.parameterEle[:,1][ind])
                distance = float(self.parameterEle[:,2][ind])
                beta1 = float(self.parameterEle[:,3][ind])
                for dis in self.dis_component[i]:
                    referenceEle -= energy*(2-(dis/distance)**beta1)
                if Math.absoluteErr(referenceEle, 0) < 0.01:
                    return 'Trained'
        if len(unknow) == 1:
            tem = 0
            for undis in dis_unknow[0]:
                tem += float(2-(undis/Rphi)**beta)
            tem_E = referenceEle/tem
            if len(self.parameterEle[0]) == 0:
                with open(self.path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([unknow[0], tem_E, Rphi, beta])
            else:
                with open(self.path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([unknow[0], tem_E, Rphi, beta])
            with open(self.path, newline='') as f:
                rows = csv.reader(f, delimiter=',')
                eleParameter = np.asarray(list(rows))
            self.parameterEle = eleParameter
            return [unknow[0], tem_E, Rphi, int(beta)]
    def run(self, Rphi=4, Rrange=1, beta=1, train_mode=False, referenceEle=None):
        if train_mode and referenceEle==None:
            return 'If train_mode is used, referenceEle must be an integer.'
        self.signEle()
        return self.energyEle(Rphi, Rrange, beta, train_mode, referenceEle)
        
        
class exchange:
    def __init__(self, frag2D, fragStruct3D, basis):
        path = './componentDatabase/exc_' + basis + '.csv'

        with open(path, newline='') as f:
            rows = csv.reader(f, delimiter=',')
            excParameter = np.asarray(list(rows))
        self.path = path
        self.frag = frag2D
        self.fragStruct = fragStruct3D
        self.parameterExc = excParameter
    def signExc(self):
        sign2d = []
        AB = ['A', 'B']
        flag = True
        for frag1d in self.frag:
            for i, f in enumerate(frag1d):
                if f not in AB:
                    startind = i
                    flag = False
                    break
            if flag:
                startind = 0
            PN = 1
            
            l, r = 0, len(frag1d)-1
            lflag, rflag = startind, startind
            sign = [PN]
            while True:
                PN *= -1
                rflag += 1
                lflag -= 1
                if rflag > r and lflag < l:
                    break
                if rflag <= r:
                    sign.append(PN)
                if lflag >= l:
                    sign.insert(0, PN)
            sign2d.append(sign)
        print(sign2d)
        print(self.frag)
        self.sign2d = sign2d
        
    def energyExc(self, Rphi, Rrange, beta, train_mode, referenceExc):
        component, count_component, dis_component = [], [], []
        for i, i_frag in enumerate(self.fragStruct[0]):
            for j, j_frag in enumerate(self.fragStruct[1]):
                if self.sign2d[0][i] * self.sign2d[1][j] > 0:
                    dis = Math.distance3D(i_frag, j_frag)
                    err = Math.absoluteErr(dis, Rphi)
                    if err < Rrange:
                        if self.frag[0][i] < self.frag[1][j]:
                            com = self.frag[0][i] + self.frag[1][j]
                        else:
                            com = self.frag[1][j] + self.frag[0][i]
                        if com not in component:
                            component.append(com)
                            count_component.append(1)
                            dis_component.append([dis])
                        else:
                            com_ind = component.index(com)
                            count_component[com_ind] += 1
                            dis_component[com_ind].append(dis)
        
        print('Exc_component: ', component, '\n', 
              'Exc_component_quantity: ', count_component, '\n',
              'Exc_component_distance: ', dis_component, '\n')
        if train_mode:
            self.component = component
            self.count_component = count_component
            self.dis_component = dis_component
            print(self.modeling(Rphi, Rrange, beta, referenceExc))
            return self.parameterExc
        
        total = 0
        for i, comp in enumerate(component):
            ind = self.parameterExc[:,0].tolist().index(comp)
            energy = float(self.parameterExc[:,1][ind])
            distance = float(self.parameterExc[:,2][ind])
            beta1 = float(self.parameterExc[:,3][ind])
            for dis in dis_component[i]:
                e = energy*(2-(dis/distance)**beta1)
                total += e
        return total
    
    def modeling(self, Rphi, Rrange, beta, referenceExc):
        try:    
            complist = self.parameterExc[:,0].tolist() 
        except:
            complist = []
        unknow, dis_unknow = [], []
        for i, comp in enumerate(self.component):
            if comp not in complist:
                unknow.append(comp)
                dis_unknow.append(self.dis_component[i])
            else:
                ind = self.parameterExc[:,0].tolist().index(comp)
                energy = float(self.parameterExc[:,1][ind])
                distance = float(self.parameterExc[:,2][ind])
                beta1 = float(self.parameterExc[:,3][ind])
                for dis in self.dis_component[i]:
                    referenceExc -= energy*(2-(dis/distance)**beta1)
                if Math.absoluteErr(referenceExc, 0) == 0:
                    return 'Trained'
        if len(unknow) == 1:
            tem = 0
            for undis in dis_unknow[0]:
                tem += float(2-(undis/Rphi)**beta)
            tem_E = referenceExc/tem
            try:
                if len(self.parameterExc[0]) == 0:
                    with open(self.path, 'w', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([unknow[0], tem_E, Rphi, beta])
                else:
                    with open(self.path, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([unknow[0], tem_E, Rphi, beta])
            except:
                with open(self.path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([unknow[0], tem_E, Rphi, beta])
            with open(self.path, newline='') as f:
                rows = csv.reader(f, delimiter=',')
                excParameter = np.asarray(list(rows))
            self.parameterExc = excParameter
            return [unknow[0], tem_E, Rphi, float(beta)]
        
    def run(self, Rphi=4, Rrange=1, beta=1, train_mode=False, referenceExc=None):
        if train_mode and referenceExc==None:
            return 'If train_mode is used, referenceExc must be an integer.'
        self.signExc()
        return self.energyExc(Rphi, Rrange, beta, train_mode, referenceExc)        
        
class induction:
    def __init__(self, frag2D, fragStruct3D, basis):
        path = './componentDatabase/ind_' + basis + '.csv'
        with open(path, newline='') as f:
            rows = csv.reader(f, delimiter=',')
            indParameter = np.asarray(list(rows))
        self.path = path
        self.frag = frag2D
        self.fragStruct = fragStruct3D
        self.parameterInd = indParameter
    
    def energyInd(self, Rphi, Rrange, train_mode, referenceInd):
        component, count_component, dis_component = [], [], []
        for i, i_frag in enumerate(self.fragStruct[0]):
            for j, j_frag in enumerate(self.fragStruct[1]):
                dis = Math.distance3D(i_frag, j_frag)
                err = Math.absoluteErr(dis, Rphi)
                if err < Rrange:
                    if self.frag[0][i] < self.frag[1][j]:
                        com = self.frag[0][i] + self.frag[1][j]
                    else:
                        com = self.frag[1][j] + self.frag[0][i]
                    if com not in component:
                        component.append(com)
                        count_component.append(1)
                        dis_component.append([dis])
                    else:
                        com_ind = component.index(com)
                        count_component[com_ind] += 1
                        dis_component[com_ind].append(dis)
        
        print('Ind_component: ', component, '\n', 
              'Ind_component_quantity: ', count_component, '\n',
              'Ind_component_distance: ', dis_component, '\n')
        if train_mode:
            self.component = component
            self.count_component = count_component
            self.dis_component = dis_component
            print(self.modeling(Rphi, Rrange, referenceInd))
            return self.parameterInd
        total = 0
        for i, comp in enumerate(component):
            ind = self.parameterInd[:,0].tolist().index(comp)
            energy = float(self.parameterInd[:,1][ind])
            distance = float(self.parameterInd[:,2][ind])
            for dis in dis_component[i]:
                e = energy*(2-(dis/distance)**3)
                total += e
        return total
    
    def modeling(self, Rphi, Rrange, referenceInd):
        try:    
            complist = self.parameterInd[:,0].tolist() 
        except:
            complist = []
        unknow, dis_unknow = [], []
        for i, comp in enumerate(self.component):
            if comp not in complist:
                unknow.append(comp)
                dis_unknow.append(self.dis_component[i])
            else:
                ind = self.parameterInd[:,0].tolist().index(comp)
                energy = float(self.parameterInd[:,1][ind])
                distance = float(self.parameterInd[:,2][ind])
                for dis in self.dis_component[i]:
                    referenceInd -= energy*(2-(dis/distance)**3)
                if Math.absoluteErr(referenceInd, 0) == 0:
                    return 'Trained'
        if len(unknow) == 1:
            tem = 0
            for undis in dis_unknow[0]:
                tem += float(2-(undis/Rphi)**3)
            tem_E = referenceInd/tem
            if len(self.parameterInd[0]) == 0:
                with open(self.path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([unknow[0], tem_E, Rphi])
            else:
                with open(self.path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([unknow[0], tem_E, Rphi])
            with open(self.path, newline='') as f:
                rows = csv.reader(f, delimiter=',')
                indParameter = np.asarray(list(rows))
            self.parameterInd = indParameter
            return [unknow[0], tem_E, Rphi]
    def run(self, Rphi=4, Rrange=1, train_mode=False, referenceInd=None):
        if train_mode and referenceInd == None:
            return 'If train_mode is used, referenceEle must be an integer.'
        return self.energyInd(Rphi, Rrange, train_mode, referenceInd)        
        
class dispersion:
    def __init__(self, frag2D, fragStruct3D, basis):
        path = './componentDatabase/disp_' + basis + '.csv'
        with open(path, newline='') as f:
            rows = csv.reader(f, delimiter=',')
            dispParameter = np.asarray(list(rows))
        self.path = path
        self.frag = frag2D
        self.fragStruct = fragStruct3D
        self.parameterDisp = dispParameter
    
    def energyDisp(self, Rrange, train_mode, referenceDisp):
        component1, component2, dis_component1, dis_component2 = [], [], [], []
        check1, check2 = [], []
        for i, i_frag in enumerate(self.fragStruct[0]):
            for j, j_frag in enumerate(self.fragStruct[1]):
                dis = Math.distance3D(i_frag, j_frag)
                if dis < Rrange:
                    if i not in check1:
                        check1.append(i)
                        component1.append(self.frag[0][i])
                        dis_component1.append(dis)
                    if j not in check2:
                        check2.append(j)
                        component2.append(self.frag[1][j])
                        dis_component2.append(dis)
        
        print(component1, '\n', component2)
        volmunDF = database.segmentVolumeDatabase()
        molecularVolume1, molecularVolume2 = 0, 0
        for c1 in component1:
            if c1 == 'D' and len(component1) == 1:
                molecularVolume1 = 4.62
                break
            elif c1 == 'F' and len(component2) == 1:
                molecularVolume1 = 3.56
                break
            vol = float(volmunDF.loc[volmunDF['segment'] == c1]['volume'])
            if molecularVolume1 == 0:
                molecularVolume1 += vol
                continue
            molecularVolume1 += np.pi*(0.9**2)*1.5 + vol
           
        for c2 in component2:
            if c2 == 'D' and len(component2) == 1:
                molecularVolume2 = 4.62
                break
            elif c2 == 'F' and len(component2) == 1:
                molecularVolume2 = 3.56
                break
            vol = float(volmunDF.loc[volmunDF['segment'] == c2]['volume'])
            if molecularVolume2 == 0:
                molecularVolume2 += vol
                continue
            molecularVolume2 += np.pi*(0.9**2)*1.5 + vol
        if train_mode:
            print(molecularVolume1,molecularVolume2)
            molvol = np.sqrt(molecularVolume1*molecularVolume2)
            return [molvol, referenceDisp], component1
        
        nonAB = ['C', 'D', 'E', 'F', 'G', 'H']
        flag1, flag2 = True, True
        for non in nonAB:
            if non in component1:
                ind1 = self.parameterDisp[:,0].tolist().index(non)
                flag1 = False
                break
        if flag1:
            ind1 = self.parameterDisp[:,0].tolist().index('AB')
        refV1 = float(self.parameterDisp[:,1][ind1])
        refE1 = float(self.parameterDisp[:,2][ind1])
        alpha1 = float(self.parameterDisp[:,3][ind1])
        for non in nonAB:
            if non in component2:
                ind2 = self.parameterDisp[:,0].tolist().index(non)
                flag2 = False
                break
        if flag2:
            ind2 = self.parameterDisp[:,0].tolist().index('AB')
        refV2 = float(self.parameterDisp[:,1][ind2])
        refE2 = float(self.parameterDisp[:,2][ind2])
        alpha2 = float(self.parameterDisp[:,3][ind2])
#         print(molecularVolume1, molecularVolume2, refE1)
#         print(refV1, refV2, refE2)
        e = -1*((molecularVolume1/refV1)**alpha1*(molecularVolume2/refV2)**alpha2)**0.5*(refE1*refE2)**0.5
        return e
    def modeling(self, volumn_refDisp_3D, component3D, maxerr):
        vr3d = np.array(volumn_refDisp_3D)
        volSet = vr3d[:,0]
        dispSet = vr3d[:,1]
        refVol = min(volSet)
        refE = dispSet[volSet.tolist().index(refVol)]
        minavg = 100
        for a in range(300):
            errlist = []
            for i, vol in enumerate(volSet):
                ratioV = (vol/refVol)**(a/100)
                predDisp = ratioV*refE
                err = abs(predDisp - dispSet[i])
                errlist.append(err)
            check = np.where(np.array(errlist) > maxerr)
            if len(check[0]) > 0:
                continue
            if sum(errlist) / len(dispSet) < minavg:
                alpha = a/100
                errResult = errlist
                minavg = sum(errlist) / len(dispSet)
        print(alpha, '\n', errResult)
        nonAB = ['C', 'D', 'E', 'F', 'G', 'H']
        ABflag = True
        for i, non in enumerate(nonAB):
            if non in component3D[-1]:
                ABflag = False
                componentName = non
                break
        if ABflag:
            componentName = 'AB'
        try:
            if len(self.parameterDisp[0]) == 0:
                with open(self.path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([componentName, refVol, refE, alpha])
            else:
                with open(self.path, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([componentName, refVol, refE, alpha])
        except:
            with open(self.path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([componentName, refVol, refE, alpha])
        with open(self.path, newline='') as f:
            rows = csv.reader(f, delimiter=',')
            dispParameter = np.asarray(list(rows))
        self.parameterDisp = dispParameter
        return [componentName, refVol, refE, alpha]
    def run(self, Rrange=5, train_mode=False, referenceDisp=None):
        return self.energyDisp(Rrange, train_mode, referenceDisp)
        
        
        
        
               
        