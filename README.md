
###A tool to estimate Instantaneous Reproduction Number(Rt) 
cal_r()provide real-time estimation of time-varying distribution of 
Rt and infected numbers from a range of epidemic observations (e.g., number of onsets and confirmed cases).

##Usage
dart = DARt(GT,D_s,Filename)

cal_r()

##Arguments
GT: generation time distribution

D_s: delay time distribution

Filename:input file

##References:
This tool is described in the following paper:

##Examples

    GT = [0, 0, 0.165720874545241, 0.226350757019051, 0.245007574714227, 0.213515210247327,
          0.149405583474155]  # 1:7; non-zero index 3:7
    D_s = np.array([0, 0, 0, 0, 0, 0, 0.0996906, 0.1130266, 0.1143032, 0.1069238, 0.0937167999999999])
    dart = DARt(GT=GT, D_s=D_s, filename='uk_report_1011.csv')
    dart.cal_r()