

#%% read the files

def main():
	with open("Data/M1.txt",'r') as f:
		M = [list(map(float,item.split())) for item in f.read().split('\n')]
	with open("Data/M2.txt",'r') as f:
		N = [ list(map(float,item.split())) for item in f.read().split('\n')]
	# print(M)
	# print(N)
	print(shape(M))




# main()
# %%





def shape(M):
	def shapeHapler(Mat):
		print(Mat)
		if(type(Mat)==list):
			l = len(Mat)
			s = shapeHapler(Mat[0])
			return  [l] + s
		else:
			return []
	return shapeHapler(M)
		
main()






# %%
