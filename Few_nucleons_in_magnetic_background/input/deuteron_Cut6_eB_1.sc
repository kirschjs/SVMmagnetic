пн фев 18 22:48:05 IST 2019
	 open input file ./input/deuteron_Cut6_eB_1.inp

   npar=2
!  masses
   xm(1)=1.  xm(2)=1.         
!  charges
!  z( 1)=0.  z( 2)=0. 
!
!   isospin configurations
    nisc=2
    cisc(1)= 1 iso(1,1)=2 iso(2,1)=1 
    cisc(2)=-1 iso(1,2)=1 iso(2,2)=2 
!
!   spin configurations
    nspc=2
    cspc(1)= 1 isp(1,1)=1 isp(2,1)=2 
    cspc(2)= 1 isp(1,2)=2 isp(2,2)=1 
!
    h2m=41.47
    irand=-1776 ibf=1
  eB=1 
    mm0=10  kk0=10  mnb=400
    bmin=0.01  bmax=6
!
! Pionless EFT LO Potential
! number of potential terms and operators
  npt=1  nop=3
 vpot(1,1)=-1036.326942     apot(1,1)=9   bpot(1,1)=0  npot(1,1)=0
 vpot(1,3)=-54.2653923      apot(1,3)=9   bpot(1,3)=0  npot(1,3)=0




	 Initialize SVM 
	 Initialize MatrixElement
	 Start SVM iters

	 iter =    1     E =    19.62078641    dE =     0.50000000 
	 iter =    2     E =    18.41048758    dE =     0.06573964 
	 iter =    3     E =    16.55938900    dE =     0.11178544 
	 iter =    4     E =    10.53091245    dE =     0.57245529 
	 iter =    5     E =    10.17953345    dE =     0.03451818 
   more eigenvalues=  39.7893  169.466  422.033  
	 iter =    6     E =     9.08461083    dE =     0.12052499 
	 iter =    7     E =     8.23389740    dE =     0.10331844 
	 iter =    8     E =     8.00537590    dE =     0.02854601 
	 iter =    9     E =     7.88883875    dE =     0.01477241 
	 iter =   10     E =     7.20132477    dE =     0.09547049 
   more eigenvalues=  28.5183  64.92  142.214  
	 iter =   11     E =     7.04850442    dE =     0.02168124 
	 iter =   12     E =     6.28495686    dE =     0.12148812 
	 iter =   13     E =     5.95767212    dE =     0.05493500 
	 iter =   14     E =     5.81002700    dE =     0.02541212 
	 iter =   15     E =     5.30605757    dE =     0.09498002 
   more eigenvalues=  28.2125  60.0846  102.363  
	 iter =   16     E =     5.15617494    dE =     0.02906857 
	 iter =   17     E =     5.10802609    dE =     0.00942612 
	 iter =   18     E =     5.08773610    dE =     0.00398802 
	 iter =   19     E =     5.03904429    dE =     0.00966291 
	 iter =   20     E =     5.00904082    dE =     0.00598986 
   more eigenvalues=  27.8577  48.4645  71.8356  
	 iter =   21     E =     4.97041261    dE =     0.00777163 
	 iter =   22     E =     4.95106283    dE =     0.00390821 
	 iter =   23     E =     4.73062378    dE =     0.04659831 
	 iter =   24     E =     4.63086198    dE =     0.02154281 
	 iter =   25     E =     4.59058213    dE =     0.00877445 
   more eigenvalues=  27.8174  48.2968  71.7858  
	 iter =   26     E =     4.40933878    dE =     0.04110443 
	 iter =   27     E =     4.37619112    dE =     0.00757455 
	 iter =   28     E =     4.34541025    dE =     0.00708353 
	 iter =   29     E =     4.31400936    dE =     0.00727882 
	 iter =   30     E =     4.27852253    dE =     0.00829418 
   more eigenvalues=  27.7333  47.714  70.0497  
	 iter =   31     E =     4.26643047    dE =     0.00283423 
	 iter =   32     E =     4.26067055    dE =     0.00135188 
	 iter =   33     E =     4.22894870    dE =     0.00750112 
	 iter =   34     E =     4.20881870    dE =     0.00478282 
	 iter =   35     E =     4.20050953    dE =     0.00197813 
   more eigenvalues=  27.2854  47.6531  67.009  
	 iter =   36     E =     4.19332806    dE =     0.00171260 
	 iter =   37     E =     4.18437577    dE =     0.00213945 
	 iter =   38     E =     4.17841667    dE =     0.00142616 
	 iter =   39     E =     4.16673794    dE =     0.00280285 
	 iter =   40     E =     4.16182329    dE =     0.00118089 
   more eigenvalues=  26.9217  45.1206  59.6271  
	 iter =   41     E =     4.15614666    dE =     0.00136584 
	 iter =   42     E =     4.14197522    dE =     0.00342142 
	 iter =   43     E =     4.13178804    dE =     0.00246556 
	 iter =   44     E =     4.07418960    dE =     0.01413740 
	 iter =   45     E =     4.04576963    dE =     0.00702462 
   more eigenvalues=  26.9103  44.9678  59.6137  
	 iter =   46     E =     4.03960985    dE =     0.00152484 
	 iter =   47     E =     4.02107092    dE =     0.00461045 
	 iter =   48     E =     4.01172285    dE =     0.00233019 
	 iter =   49     E =     3.99579843    dE =     0.00398529 
	 iter =   50     E =     3.99191859    dE =     0.00097192 
   more eigenvalues=  26.9007  44.9502  59.6133  
	 iter =   51     E =     3.97875463    dE =     0.00330856 
	 iter =   52     E =     3.97279964    dE =     0.00149894 
	 iter =   53     E =     3.96815211    dE =     0.00117121 
	 iter =   54     E =     3.96091273    dE =     0.00182771 
	 iter =   55     E =     3.94470774    dE =     0.00410803 
   more eigenvalues=  26.5275  43.512  58.8794  
	 iter =   56     E =     3.93598672    dE =     0.00221571 
	 iter =   57     E =     3.92434808    dE =     0.00296575 
	 iter =   58     E =     3.91903389    dE =     0.00135600 
	 iter =   59     E =     3.91100091    dE =     0.00205395 
	 iter =   60     E =     3.74045379    dE =     0.04559530 
   more eigenvalues=  26.5125  43.4371  58.8647  
	 iter =   61     E =     3.66379701    dE =     0.02092277 
	 iter =   62     E =     3.65522562    dE =     0.00234497 
	 iter =   63     E =     3.64765751    dE =     0.00207479 
	 iter =   64     E =     3.63478153    dE =     0.00354244 
	 iter =   65     E =     3.63133914    dE =     0.00094797 
   more eigenvalues=  26.503  43.3127  58.8281  
	 iter =   66     E =     3.62994799    dE =     0.00038324 
	 iter =   67     E =     3.62803737    dE =     0.00052663 
	 iter =   68     E =     3.62295705    dE =     0.00140226 
	 iter =   69     E =     3.61493944    dE =     0.00221791 
	 iter =   70     E =     3.57743494    dE =     0.01048363 
   more eigenvalues=  26.0921  39.9992  56.9916  
	 iter =   71     E =     3.26044006    dE =     0.09722457 
	 iter =   72     E =     3.25950507    dE =     0.00028685 
	 iter =   73     E =     3.25772655    dE =     0.00054594 
	 iter =   74     E =     3.25548598    dE =     0.00068825 
	 iter =   75     E =     3.25391294    dE =     0.00048343 
   more eigenvalues=  25.0646  39.0151  55.6707  
	 iter =   76     E =     3.25225695    dE =     0.00050918 
	 iter =   77     E =     3.23327456    dE =     0.00587095 
	 iter =   78     E =     3.22887292    dE =     0.00136321 
	 iter =   79     E =     3.22833099    dE =     0.00016787 
	 iter =   80     E =     3.22798464    dE =     0.00010729 
   more eigenvalues=  25.0609  39.0104  55.6625  
	 iter =   81     E =     3.22666945    dE =     0.00040760 
	 iter =   82     E =     3.22652679    dE =     0.00004421 
	 iter =   83     E =     3.22329418    dE =     0.00100289 
	 iter =   84     E =     3.22299824    dE =     0.00009182 
	 iter =   85     E =     3.22260028    dE =     0.00012349 
   more eigenvalues=  24.9836  38.9372  55.6308  
	 iter =   86     E =     3.22216957    dE =     0.00013367 
	 iter =   87     E =     3.22189055    dE =     0.00008660 
	 iter =   88     E =     3.22147754    dE =     0.00012821 
	 iter =   89     E =     3.21841994    dE =     0.00095003 
	 iter =   90     E =     3.21627433    dE =     0.00066711 
   more eigenvalues=  24.9828  38.9305  55.6217  
	 iter =   91     E =     3.21605961    dE =     0.00006676 
	 iter =   92     E =     3.21580942    dE =     0.00007780 
	 iter =   93     E =     3.21542637    dE =     0.00011913 
	 iter =   94     E =     3.21499358    dE =     0.00013462 
	 iter =   95     E =     3.21489302    dE =     0.00003128 
   more eigenvalues=  24.8354  38.0999  55.2573  
	 iter =   96     E =     3.21412238    dE =     0.00023977 
	 iter =   97     E =     3.21401718    dE =     0.00003273 
	 iter =   98     E =     3.21338009    dE =     0.00019826 
	 iter =   99     E =     3.21283193    dE =     0.00017062 
	 iter =  100     E =     3.21224873    dE =     0.00018155 
   more eigenvalues=  24.8263  38.0474  55.1959  
	 iter =  101     E =     3.21215585    dE =     0.00002892 
	 iter =  102     E =     3.21200853    dE =     0.00004587 
finding new state with lower energy failed

time=  39.34
done!
пн фев 18 22:55:58 IST 2019
rm: cannot remove ‘sbatch_script.deuteron_Cut6_eB_1.16817’: No such file or directory
