import numpy as np
import pandas as pd
from database import bondLengthDatabase

class Math:
    def deviation(measure, real):
        '''to do: 計算誤差'''
        return abs(float(measure - real) / real * 100)

    def absoluteErr(measure, real):
        '''to do: 絕對誤差(kcal/mol)'''
        return abs(float(measure - real))
    
    def distance3D(point1, point2):
        '''to do: 空間中兩點之間距離'''
        return np.sqrt((float(point1[0]) - float(point2[0]))**2 + (float(point1[1]) - float(point2[1]))**2 + (float(point1[2]) -float(point2[2]))**2)

    def midpoint(point1, point2):
        '''to do: 找中點'''
        return np.array([round((point1[0] + point2[0])/2, 3), round((point1[1] + point2[1])/2, 3), round((point1[2] + point2[2])/2, 3)])

    def simultaneous(A, B):
        '''to do: 解聯立方程 x = (A^-1)dotB'''
        l = len(B)
        B = np.array(B).reshape(l, 1)
        inv_A = np.linalg.inv(A)
        ans = inv_A.dot(B)
        return ans.reshape(l)
 
    def centroid(struct2D):
        '''to do: 找質心'''
        massDict = {'C' : 12, 'O' : 16, 'N' :  14, 'H' : 1}
        centroid2D = []

        cent = np.array([0.0, 0.0, 0.0])
        totMass = 0
        for atom in struct2D:
            m = massDict[atom[0]]
            cent += np.array([float(atom[1])*m, float(atom[2])*m, float(atom[3])*m])
            totMass += m
        cent = cent/totMass
        cent = cent.tolist()

        return cent
                
class Convert:    
    def fragmentConvert(struct2D):
        '''to do: 把原本的分子轉換成粗粒化大分子'''
        bondDF = bondLengthDatabase()
        bondName = bondDF['bond'].to_list()
        bondDis = bondDF['distance'].to_list()
        frag, fragStruct = [], []
        nextflag = False
        for i, struct in enumerate(struct2D):
            if struct[0] != 'C' or nextflag:
                nextflag = False
                continue
#             print(struct2D[i][0], struct2D[i+1][0], struct2D[i+2][0])
            if struct2D[i+1][0] != 'C' and struct2D[i+2][0] != 'C' and struct2D[i+2][0] != 'H':
                if struct2D[i+1][0] == 'N' or struct2D[i+2][0] == 'N':
                    frag.append('H')
#                     print(struct2D[i:i+3])
                    fragStruct.append(Math.centroid(struct2D[i:i+3]))
                else:
                    frag.append('G')
#                     print(struct2D[i:i+3])
                    fragStruct.append(Math.centroid(struct2D[i:i+3]))
            elif struct2D[i+1][0] == 'O' :
                dis = Math.distance3D(struct2D[i][1:], struct2D[i+1][1:])
                for dind, d in enumerate(bondDis):
                    err = Math.deviation(dis, d)
                    if err < 5:
                        BN = struct2D[i][0] + struct2D[i+1][0]
                        if BN == bondName[dind][:2]:
                            if bondName[dind][2:] == '1':
                                frag.append('E')
#                                 print(struct2D[i:i+2])
                                fragStruct.append(Math.centroid(struct2D[i:i+2]))
                            elif bondName[dind][2:] == '2':
                                frag.append('F')
#                                 print(struct2D[i:i+2])
                                fragStruct.append(Math.centroid(struct2D[i:i+2]))
            elif struct2D[i+1][0] != 'H':
                dis = Math.distance3D(struct2D[i][1:], struct2D[i+1][1:])
                for dind, d in enumerate(bondDis):
                    err = Math.deviation(dis, d)
                    if err < 5:
                        BN = struct2D[i][0] + struct2D[i+1][0]
                        if BN == bondName[dind][:2]:
                            if bondName[dind][2:] == '1':
                                if i == 0:
                                    frag.append('A')
#                                     print(struct)
                                    fragStruct.append(Math.centroid([struct]))
                                else:
                                    frag.append('B')
#                                     print(struct)
                                    fragStruct.append(Math.centroid([struct]))
                            elif bondName[dind][2:] == '2':
                                frag.append('C')
                                nextflag = True
#                                 print(struct2D[i:i+2])
                                fragStruct.append(Math.centroid(struct2D[i:i+2]))
                            elif bondName[dind][2:] == '3':
                                frag.append('D')
                                nextflag = True
#                                 print(struct2D[i:i+2])
                                fragStruct.append(Math.centroid(struct2D[i:i+2]))
            else:
                frag.append('A')

                fragStruct.append(Math.centroid([struct]))

        return frag, fragStruct
    
    
    
class Sort:
    def fragmentSort(struct3D):
        '''to do: 排序原子'''
        sortStruct = []
        bondDF = bondLengthDatabase()
        bondName = bondDF['bond'].to_list()
        bondDis = bondDF['distance'].to_list()
        for a in range(len(struct3D)):
            for st in range(len(struct3D[a])):
                if struct3D[a][st][0] == 'C':
                    Cstruct, Cind = [], []
                    nonCHstruct, nonCHind = [], []
                    for stru in struct3D[a]:
                        if stru[0] == 'C':
                            Cstruct.append(stru)
                            Cind.append(struct3D[a].index(stru))
                        elif stru[0] != 'C' and stru[0] != 'H':
                            nonCHstruct.append(stru)
                            nonCHind.append(struct3D[a].index(stru))
                    def recursiveIndex(index, Cind, nonCHind, struct2D, result = [], count = 0):
                        if count == 0:
                            result.append(struct2D[index])
                            if len(nonCHind) > 0:
                                for non in nonCHind:
#                                     print(struct2D[index][0], struct2D[non][0])
                                    nondis = Math.distance3D(struct2D[index][1:], struct2D[non][1:])
#                                     print(nondis, '-----------------')
                                    for dind, d in enumerate(bondDis):
                                        err = Math.deviation(nondis, d)
                                        if err < 5:
                                            BN = struct2D[index][0] + struct2D[non][0]
                                            if BN == bondName[dind][:2]:
                                                result.append(struct2D[non])
                            Cind.remove(index)
                        count += 1
                        if len(Cind) > 0:
                            for i, s in enumerate(struct2D):
                                if s[0] == 'C':
                                    if i in Cind:
                                        dis = Math.distance3D(struct2D[index][1:], s[1:])
                                        err1 = Math.deviation(dis, 1.52)
                                        err2 = Math.deviation(dis, 1.33)
                                        err3 = Math.deviation(dis, 1.21)
                                        if err1 < 5 or err2 < 5 or err3 < 5:
                                            result.append(s)
                                            if len(nonCHind) > 0:
                                                for non in nonCHind:
#                                                     print(struct2D[index][0], struct2D[non][0])
                                                    nondis = Math.distance3D(s[1:], struct2D[non][1:])
#                                                     print(nondis, '-----------------')
                                                    for dind, d in enumerate(bondDis):
                                                        err = Math.deviation(nondis, d)
#                                                         print(err)
                                                        if err < 5:
                                                            BN = s[0] + struct2D[non][0]
                                                            if BN == bondName[dind][:2]:
                                                                result.append(struct2D[non])
                                            Cind.remove(i)
                                            return recursiveIndex(i, Cind, nonCHind, struct2D, result, count)
                            return False
                        else: return result
                    res = recursiveIndex(st, Cind, nonCHind, struct3D[a])
                    if res:break
            try:
                for stru in struct3D[a]:
                    if stru[0] == 'H':
                        res.append(stru)
                sortStruct.append(res)
            except:
                return struct3D
        return sortStruct
    
    def dicectionSort(frag2D, fragStruct3D):
        '''to do: 方向比對'''
        vector = []
        for fst in fragStruct3D:
            vec = np.array([0,0,0], dtype='float64')
            for i in range(len(fst)-1):
                v = np.array(fst[i+1]) - np.array(fst[i])
                vec += v
            vector.append(vec)
        if np.dot(vector[0], vector[1]) >= 0:
            return frag2D, fragStruct3D
        else:
            fragStruct3D[1] = np.flip(fragStruct3D[1], 0).tolist()
            frag2D[1] = np.flip(frag2D[1], 0).tolist()
            return frag2D, fragStruct3D

    
    
    
    
    

        
        
        
        
        
        
        
        
        
        
    
    
    
    