import os
import xlwt 

def isfloat(x):
    try:
        float(x)
        return True
    except:
        return False

book=xlwt.Workbook(encoding='utf-8')
sheet=book.add_sheet('sheet1')

filepath=r'D:\User\Spectra\GEM20\\'
list_file=os.listdir(filepath)
row=0
for file in list_file:
    print(file[-4:])
    if file[-4:] == '.Spe':
        inp=file
        line=0
        file_open=open(filepath+inp,'r',encoding='unicode_escape')
        file_read=file_open.read()
        ind=file_read.index('$DATA:')
        ind+=14
        # print(file_read[ind:ind+100])
        sheet.write(line,row,inp)
        line=1
        while True:
            cts=file_read[ind:ind+9]
            isf=isfloat(cts)
            if isf == True:
                cts=float(cts)
                sheet.write(line,row,cts)
                ind+=9
                line+=1
            elif '$ROI:' in cts:
                break
            else:
                continue
        row+=1
book.save(r'D:\User\Spectra\GEM20\total.xls')
