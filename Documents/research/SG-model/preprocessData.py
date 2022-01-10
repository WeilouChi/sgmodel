import numpy as np
import os
import sys
import shutil
from calculate import Math

def start_stop(data):
    '''to do: 找出開始和結束的index'''
    ind = '0'
    for i in range(len(data)):
        '''to do: 去找input file在第幾行'''
        
        if data[i][0:20] == '  ==> Input File <==' :
            ind = i
            break
    if ind == '0':
        print('INPUT FILE ERROR: Sapt output file has no input file.')
        sys.exit()

    start_ind, stop_ind = '0', '0'        
    flag = False
    for i in range(ind+3, len(data)):
        '''to do: 去找結構開始和結束的index'''
        start = data[i].split()
        try:
            if start[0] == 'molecule' and start[-1] == '{':
                start_ind = i
                flag = True
            elif flag:
                if start[-1] == '}':
                    stop_ind = i
                    break
        except:
            continue
    if start_ind == '0' or stop_ind == '0':
        print('INPUT FILE ERROR: Sapt output file format does not match.')
        sys.exit()
    return start_ind, stop_ind  


def find_energy(data):
    '''to do : find Sapt0 result'''

    ind2 = '0'
    for i in data:
        if i[:16] == '    SAPT Results':
            ind2 = data.index(i)
            break
    if ind2 == '0':
        print('INPUT FILE ERROR: Sapt output file has no Sapt result.')
        sys.exit()

    energy = np.empty([0,4])
    flag = False
    for i in range(ind2, len(data)):
        '''to do: 找出四種力和總和'''
        ere = data[i].split()
        if flag:
            break
        try:
            if ere[0] == 'Electrostatics':
                energy = np.append(energy, float(ere[3]))
            elif ere[0] == 'Exchange':
                energy = np.append(energy, float(ere[3]))
            elif ere[0] == 'Induction':
                energy = np.append(energy, float(ere[3]))
            elif ere[0] == 'Dispersion':
                energy = np.append(energy, float(ere[3]))
            elif ere[0] == 'Total' and ere[1] == 'SAPT0':
                energy = np.append(energy, float(ere[4]))
                flag = True
        except:
            continue
    return energy


def sort_molecularWeight(molecularWeight, spec):
    '''to do : sort 分子量大到小'''
    tem = []
    result = []
    for i in spec:
        tem.append(molecularWeight[i])
    tem = sorted(tem, reverse=True)
    for i in tem:
        for j in molecularWeight:
            if i == molecularWeight[j]:
                result.append(j)
    return result


def find_structure(start_ind, stop_ind, data):
    '''to do: 找出資料中的結構'''
    molecularWeight = {'C': 12, 'O' : 16, 'N' : 14, 'H' : 1}
    atoms = [0] * len(molecularWeight)
    spec = [i for i in molecularWeight]
    spec = sort_molecularWeight(molecularWeight, spec)
    all_molecular = []
    all_structure = []
    structure = []
    for i in range(start_ind, stop_ind):
        '''to do: 分析這個結構是dimer or bimer, 把原子存進atoms, 並按照分子量由大到小, 再把分子結構取出'''
        atom = data[i].split()
        try:
            if atom[0] == '--':
                all_molecular.append(atoms)
                all_structure.append(structure)
                structure = []
                atoms = [0] * len(molecularWeight)
            if atom[0] in molecularWeight:
                structure.append(data[i][:len(data[i])-2].split())
                spec_index = spec.index(atom[0])
                atoms[spec_index] += 1
                
        except:
            continue
    all_molecular.append(atoms)
    all_structure.append(structure)

    return all_molecular, all_structure, spec


def sort_structure(structure):
    '''to do: 把相同的原子排在一起'''
    result = []
    spec_id = ['C', 'O', 'N', 'H']
    for Id in spec_id:
        for i in structure:
            if i[0] == Id:
                result.append(i)
    return result

def getfinal_stru(all_structure):
    '''to do: 找出最後要使用的data'''
    sortedStruc = []
    for i in all_structure:
        '''to do: 排序所有的構型'''
        tem = sort_structure(i)
        sortedStruc.append(tem)
    return sortedStruc

def judgeMolecular(all_molecular, spec):
    '''to do: 判斷是什麼構型'''
    if all_molecular[0] == all_molecular[1]:
        '''to do : 判斷是Dimer or Bimer'''
        molecularType = 'Dimer'
    else:
        molecularType = 'Bimer'
        
    orgCom = []

    for sin in all_molecular:
        if sin[spec.index('C')] * 2 + 2 == sin[spec.index('H')] and sin[spec.index('O')] == 0 and sin[spec.index('N')] == 0:
            orgCom.append('Alkane')
                                                                
        elif sin[spec.index('C')] * 2 == sin[spec.index('H')] and sin[spec.index('O')] == 0 and sin[spec.index('N')] == 0:
            orgCom.append('Alkene')
                                                                    
        elif sin[spec.index('C')] * 2 - 2 == sin[spec.index('H')] and sin[spec.index('O')] == 0 and sin[spec.index('N')] == 0:
            orgCom.append('Alkyne')
                                                                        
        elif sin[spec.index('C')] * 2 + 2 == sin[spec.index('H')] and sin[spec.index('O')] == 1 and sin[spec.index('N')] == 0:
            orgCom.append('Alcohol')
                                                                        
        elif sin[spec.index('C')] * 2 == sin[spec.index('H')] and sin[spec.index('O')] == 1 and sin[spec.index('N')]== 0:
            orgCom.append('Aldehyde')                                                                
                                
        elif sin[spec.index('C')] * 2 + 4 == sin[spec.index('H')] and sin[spec.index('O')] == 1 and sin[spec.index('N')] == 0:
            orgCom.append('Ketone')       
                                                                        
        elif sin[spec.index('C')] * 2 == sin[spec.index('H')] and sin[spec.index('O')] == 2 and sin[spec.index('N')] == 0:
            orgCom.append('Acid')
        
        elif sin[spec.index('C')] * 2 + 1 == sin[spec.index('H')] and sin[spec.index('O')] == 1 and sin[spec.index('N')] == 1:
            orgCom.append('Amide')
    
    return molecularType, orgCom

def find_basis(data):
    '''to do: 找出需要的基底 只能有jdz or jtz'''
    for i in range(len(data)):
        bas = data[i].split()
        #print(bas)
        try:
            if bas[0] == 'basis':
                if bas[1].lower() == 'jun-cc-pvdz' or bas[1].lower() == 'jdz':
                    return 'jdz'
                elif bas[1].lower() == 'jun-cc-pvtz' or bas[1].lower() == 'jtz':
                    return 'jtz'
                else:
                    print('INPUT FILE ERROR: BASIS DO NOT MATCH.')
                    sys.exit()
        except:
            continue

     
    
def dataMove():
    dataPath = './initialize_data'
    fileNamelist = os.listdir(dataPath)
    nowPath = os.getcwd()
    for file in fileNamelist:
        beforePath = os.path.join(nowPath, "initialize_data")
        before = os.path.join("initialize_data", file)
        pdgod = getOutputData()
        orgType = pdgod.dataOrgcom(file)
        basis = pdgod.dataBasis(file)
        try:
            org = orgType[0] + '-' + orgType[1]
            afterPath = os.path.join(nowPath, "molecularStructure")
            afterPath = os.path.join(afterPath, basis)
            afterPath = os.path.join(afterPath, org)
            after = os.path.join(afterPath, file)
            shutil.copyfile(before, after)
        except FileNotFoundError:
            
            org = orgType[1] + '-' + orgType[0]
            afterPath = os.path.join(nowPath, "molecularStructure")
            afterPath = os.path.join(afterPath, basis)
            afterPath = os.path.join(afterPath, org)
            after = os.path.join(afterPath, file)
            shutil.copyfile(before, after)
        os.remove(before)


def dataName(path='./initialize_data'):
    file_name = []
    fileNamelist = os.listdir(path)
    for file in fileNamelist:
        file_name.append(file)
    return file_name
        
class getOutputData:
    def __init__(self, path='./initialize_data'):
        self.path = path
        
    def dataEnergy(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        start_ind, stop_ind = start_stop(data)
        energy = find_energy(data)
        energy = energy.tolist()
        return energy

    def dataStructure(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        start_ind, stop_ind = start_stop(data)
        all_molecular, all_structure, spec = find_structure(start_ind, stop_ind, data)
        sortedStruc = getfinal_stru(all_structure)
        return sortedStruc

    def dataMolecular(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        start_ind, stop_ind = start_stop(data)
        all_molecular, all_structure, spec = find_structure(start_ind, stop_ind, data)
        return all_molecular
    
    def dataSpec(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        start_ind, stop_ind = start_stop(data)
        all_molecular, all_structure, spec = find_structure(start_ind, stop_ind, data)
        return spec

    def dataMoltype(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
            start_ind, stop_ind = start_stop(data)
            all_molecular, all_structure, spec = find_structure(start_ind, stop_ind, data)
            molecularType, orgCom = judgeMolecular(all_molecular, spec)
        return molecularType

    def dataOrgcom(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
            start_ind, stop_ind = start_stop(data)
            all_molecular, all_structure, spec = find_structure(start_ind, stop_ind, data)
            molecularType, orgCom = judgeMolecular(all_molecular, spec)
        return orgCom

    def dataBasis(self, file):
        with open (self.path + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        basis = find_basis(data)
        return basis



class getIntputData:
    def dataStructure(file):
        dataPath = './initialize_data'
        start_ind = []
        stop_ind = []
        final_structure = []
        with open (dataPath + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        for i in range(len(data)):
            dat = data[i].split()
            try:
                if dat[0] == '0' and dat[1] == '1' and len(dat) == 2:
                    start_ind.append(i)
                elif dat[0] == '--':
                    stop_ind.append(i)
            except:
                pass
        sturc1 = []
        for i in range(start_ind[0]+1, stop_ind[0]):
            dat = data[i].split()
            try:
                if dat[0]:
                    sturc1.append(dat)
            except:
                pass
        sturc2 = []
        for i in range(start_ind[1]+1, len(data)):
            dat = data[i].split()
            try:
                if dat[0] == '}':
                    break
                elif dat[0] == 'units' or dat[0] == 'no':
                    break
            except:
                 dat.append(None)
            if dat[0] != None:
                sturc2.append(dat)
        final_structure = [sturc1, sturc2]
        return final_structure
    
    def dataBasis(file):
        dataPath = './initialize_data'
        with open (dataPath + '/' + file , 'r' ) as init_data:
            data = init_data.readlines()
        basis = find_basis(data)
        return basis

class writeData:
    def gaussianView(name, struct3D):
        '''to do: make .gjf file.'''
        gaussian = ["%chk=C:\\Users\\User\\Desktop\\gaussian\\pentane\\pentane.chk\n",'# hf/3-21g geom=connectivity\n','\n',
            'Title Card Required\n','\n','0 1\n']
        for i in struct3D:
            for j in i:
                row = j[0] + '         ' + j[1] + '   ' + j[2] + '   ' + j[3] + '\n'
                gaussian.append(row)
        with open('.\\gaussianView\\%s.gjf' % name, 'w') as gasdata:
            gasdata.writelines(gaussian)

    def gViewToSaptInput(gViewfolderPath, saptfolderPath):
        '''to do: .gjf to sapt0 input file.'''
        import os
        for file in os.listdir(gViewfolderPath):
            inp = ["# Any line starting with the # character is a comment line\n", '#! Sample  computation\n', '\n', 'memory 22 gb\n', '\n', 'molecule Dimer {\n','\n','0 1\n']
            path = os.path.join(gViewfolderPath, file)
            if file[-3:] == 'gjf':
                with open (path, 'r') as intdata:
                    data = intdata.readlines()
            else:continue
            flag, flag2 = False, False
            outdata = []
            outdata2 = []
            change = [None]
            for d in data:
                dat = d.split()
                try:
                    if dat[0] == '0' and dat[1] == '1':
                        flag = True
                        continue
                except:pass
                if flag and flag2==False:
                    try:
                        if dat[0] and dat[0] != '?s':
                            if change[0] == 'H' and dat[0] == 'C':
                                flag2 = True
                                outdata2.append(d)
                            else:
                                change[0] = dat[0]
                                outdata.append(d)
                    except:pass
                elif flag and flag2:
                    try:
                        if dat[0] and dat[0] != '?s':
                            outdata2.append(d)
                    except:
                        break
            for o1 in outdata:
                inp.append(o1)
            inp.append('--\n')
            inp.append('0 1\n')
            for o2 in outdata2:
                inp.append(o2)
            inp.append('\n')
            inp.append('units angstrom')
            inp.append('\n')
            inp.append('no_reorient\n')
            inp.append('symmetry c1\n')
            inp.append('}\n')
            inp.append('\n')
            inp.append('set globals{\n')
            inp.append('  basis jun-cc-pvtz\n')
            inp.append('}\n')
            inp.append('\n')
            inp.append("energy('sapt0')\n")
            outpath = os.path.join(saptfolderPath, file[:-4])
            with open(outpath + '.dat' , 'w') as outdata:
                outdata.writelines(inp)
                outdata.close()
        return 'Done!'
    def SaptOutputToCsv(SaptPath='./SaptOutput', CsvPath='./molecularStructure'):
        '''to do: .out to .csv'''
        for fileName in os.listdir(SaptPath):
            print(fileName[:-3])
            with open (os.path.join(SaptPath, fileName), 'r') as d:
                data = d.readlines()
                structure0, structure1 = [],[]
                flag = True
                for ind in range(data.index('molecule Dimer {\n')+1, len(data)):
                    dat = data[ind].split()
                    try:
                        if dat[0] == '--':
                            flag = False
                    except:continue
                    try:
                        if len(dat) == 4 and float(dat[1]) and float(dat[2]) and float(dat[3]):
                            if flag:
                                structure0.append(dat)
                            else:
                                structure1.append(dat)
                    except:continue
                    if dat[0] == 'basis':
                        if dat[1].lower() == 'jun-cc-pvtz':
                            basis = 'jtz'
                        elif dat[1].lower() == 'jun-cc-pvdz':
                            basis = 'jdz'
                        print(basis)
                        break

            outfileName = fileName[:-3] + 'csv'            
            with open (os.path.join(CsvPath, outfileName), 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['molecular0'])
                writer.writerows(structure0)
                writer.writerow(['molecular1'])
                writer.writerows(structure1)


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        