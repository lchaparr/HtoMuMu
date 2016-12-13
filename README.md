# HtoMuMu

##Gen Level
Studies over the kinematics of the HtoMuMu samples.

	export SCRAM_ARCH=slc6_amd64_gcc481
	cmsrel CMSSW_7_1_20
	cd CMSSW_7_1_20/src
	cmsenv
	git cms-addpkg GeneratorInterface/Pythia6Interface
	
To convert the LHE file into Root file and Hadronize:
	mkdir Configuration/GenProduction/python
	cd Configuration/GenProduction/python
	
copy the file: test_PowhegPythiaH190_cff.py and modifies the file according your necesities.	
	scram b -j 10
	
	cd ../../..
	cmsRun Configuration/GenProduction/python/test_PowhegPythiaH190_cff.py
	
After that step you obtain a py. file test_PowhegPythiaH190_cff_py_LHE_GEN.py. Modify the path of the LHE file location.
  cmsRun test_PowhegPythiaH190_cff_py_LHE_GEN.py
  
At this point, you get a root file.

To Run the GenLevel analisis code:

	git clone https://github.com/lchaparr/HtoMuMu
	mv HtoMuMu/BasicTester.cc GeneratorInterface/Pythia6Interface/test
	scram b
	cmsRun BasicTester_cfg.py

##RecoLevel
Based on https://github.com/florez/NtupleAnalyzer/tree/master/BSM3G_TNT_simpleTemplete_Analyzer

	git clone https://github.com/florez/NtupleAnalyzer
	cd BSM3G_TNT_simpleTemplete_Analyzer
Replace *.cc and *.h

	make BSM_Analysis
	./BSM_Analysis INPUT_FILE OUTPUT_FILE
It is mandatory have an n-tuple made using the TNT Maker.

##FSR Photons
Photons from PF paticle collection. The Selector, save the relevant parameters as phi, eta, pt, energy, RelISo, etc., in the n-tuples. The code for ntuplizer (https://github.com/florez/NtupleMaker_740/tree/for_CMSSW_8X) must be modify in order to obtain the FSR photons using this selector.
