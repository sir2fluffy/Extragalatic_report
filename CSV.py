import csv
import numpy as np


class CSV:
    "convineince module I wrote to handle csvs, slow but simple"
    name = ''
    path = '.csv'
    row_count = 0
    column_count = 0
    
    row_count_data = -1
    
    
    def __init__(self,name,skip_first = True):
        path = f"{name}.csv"
        self.name = name
        self.path = path
        with open(self.path, 'r', newline='') as outfile:
            reader = csv.reader(outfile)
            row_count = sum(1 for row in reader)
            self.row_count = row_count
            self.row_count_data = row_count - 1
        with open(self.path, 'r', newline='') as outfile:
            reader = csv.reader(outfile)
            for index,row in enumerate(reader):

                self.column_count = len(row)
                return

        
    def add_row(self,*row):
        with open(self.path, 'a', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(row)
            
    def add_rows(self,rows):
        with open(self.path, 'a', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(rows)
    
    def remove_row(self,row_index):
        current = self.read_all()
        with open(self.path, 'w', newline='') as outfile:
            writer = csv.writer(outfile)     
            for index,row in enumerate(current):
                if index != row_index:
                    writer.writerow(row)
    
    def remove_rows(self,row_indexs):
        current = self.read_all()
        with open(self.path, 'w', newline='') as outfile:
            writer = csv.writer(outfile)     
            for index,row in enumerate(current):
                if not index in row_indexs:
                    writer.writerow(row)
        
    def read_all(self,skip_first = True):
        
        rows = self.row_count - int(skip_first)
        temp = np.zeros((rows,self.column_count))
        
        
        with open(self.path, 'r', newline='') as outfile:
            reader = csv.reader(outfile)
            for index,row in enumerate(reader):

                if index >= int(skip_first):
                    temp[index- int(skip_first),:] = row
                        

            return temp
    
    def read_row(self,row_index):
        with open(self.path, 'r', newline='') as outfile:
            reader = csv.reader(outfile)
            return reader[row_index]
        
    def show_info(self):
        print(f'columns: {self.column_count}')
        print(f'rows: {self.row_count}')
        print(f'file name: {self.name}')
        print(f'file path: {self.path}')




