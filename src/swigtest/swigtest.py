import os

if __name__ == '__main__':
	base = os.path.dirname(os.path.realpath(__file__))
	files = os.listdir(base)
	tests = [f for f in files if '.py' in f and 'swigtest.py' != f and 'test_' in f]
	for test in tests:
		print(14*'*****'+'\nRunning: '+test+'\n'+14*'*****'+'\n')
		os.system('python '+base+'/'+test)
