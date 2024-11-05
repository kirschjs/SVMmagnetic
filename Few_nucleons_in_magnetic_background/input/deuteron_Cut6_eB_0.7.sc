пн фев 18 22:48:05 IST 2019
	 open input file ./input/deuteron_Cut6_eB_0.7.inp

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
  eB=0.7 
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

	 iter =    1     E =    14.07013579    dE =     0.50000000 
	 iter =    2     E =    13.36912544    dE =     0.05243502 
	 iter =    3     E =    12.02399543    dE =     0.11187047 
	 iter =    4     E =     7.37799082    dE =     0.62971136 
	 iter =    5     E =     7.08663757    dE =     0.04111305 
   more eigenvalues=  33.308  168.973  416.694  
	 iter =    6     E =     6.20421208    dE =     0.14223006 
	 iter =    7     E =     5.54230840    dE =     0.11942744 
	 iter =    8     E =     5.30451256    dE =     0.04482897 
	 iter =    9     E =     5.21098774    dE =     0.01794762 
	 iter =   10     E =     4.63484517    dE =     0.12430676 
   more eigenvalues=  22.5238  55.8578  142.682  
	 iter =   11     E =     4.53995018    dE =     0.02090221 
	 iter =   12     E =     3.89111040    dE =     0.16674926 
	 iter =   13     E =     3.53902570    dE =     0.09948634 
	 iter =   14     E =     3.41212656    dE =     0.03719063 
	 iter =   15     E =     2.97126413    dE =     0.14837537 
   more eigenvalues=  22.0489  52.1076  101.081  
	 iter =   16     E =     2.84035297    dE =     0.04608975 
	 iter =   17     E =     2.77745519    dE =     0.02264583 
	 iter =   18     E =     2.75790794    dE =     0.00708771 
	 iter =   19     E =     2.73164078    dE =     0.00961589 
	 iter =   20     E =     2.70466327    dE =     0.00997444 
   more eigenvalues=  21.9244  49.4675  87.427  
	 iter =   21     E =     2.65187580    dE =     0.01990571 
	 iter =   22     E =     2.63233145    dE =     0.00742473 
	 iter =   23     E =     2.43172654    dE =     0.08249484 
	 iter =   24     E =     2.41353245    dE =     0.00753837 
	 iter =   25     E =     2.39500972    dE =     0.00773389 
   more eigenvalues=  21.6706  35.114  62.6163  
	 iter =   26     E =     2.27836425    dE =     0.05119702 
	 iter =   27     E =     2.20842590    dE =     0.03166887 
	 iter =   28     E =     2.18394206    dE =     0.01121085 
	 iter =   29     E =     2.15367323    dE =     0.01405451 
	 iter =   30     E =     2.14174636    dE =     0.00556876 
   more eigenvalues=  21.57  35.0384  62.4881  
	 iter =   31     E =     2.12948245    dE =     0.00575910 
	 iter =   32     E =     2.11159983    dE =     0.00846875 
	 iter =   33     E =     2.09515173    dE =     0.00785055 
	 iter =   34     E =     2.08741296    dE =     0.00370735 
	 iter =   35     E =     2.07494052    dE =     0.00601099 
   more eigenvalues=  21.3539  34.7714  53.7327  
	 iter =   36     E =     2.06840406    dE =     0.00316015 
	 iter =   37     E =     2.05640804    dE =     0.00583348 
	 iter =   38     E =     2.04008501    dE =     0.00800115 
	 iter =   39     E =     2.03345827    dE =     0.00325885 
	 iter =   40     E =     2.02918440    dE =     0.00210620 
   more eigenvalues=  21.3458  34.6868  53.6137  
	 iter =   41     E =     2.02128899    dE =     0.00390612 
	 iter =   42     E =     2.00645920    dE =     0.00739103 
	 iter =   43     E =     2.00143095    dE =     0.00251233 
	 iter =   44     E =     1.96149734    dE =     0.02035874 
	 iter =   45     E =     1.94055215    dE =     0.01079342 
   more eigenvalues=  21.203  34.5591  53.3975  
	 iter =   46     E =     1.93786803    dE =     0.00138509 
	 iter =   47     E =     1.89051596    dE =     0.02504716 
	 iter =   48     E =     1.87498096    dE =     0.00828542 
	 iter =   49     E =     1.86944969    dE =     0.00295877 
	 iter =   50     E =     1.86323151    dE =     0.00333731 
   more eigenvalues=  21.1646  34.5437  53.2924  
	 iter =   51     E =     1.85254545    dE =     0.00576831 
	 iter =   52     E =     1.84860518    dE =     0.00213149 
	 iter =   53     E =     1.84490629    dE =     0.00200492 
	 iter =   54     E =     1.81786274    dE =     0.01487656 
	 iter =   55     E =     1.81090743    dE =     0.00384078 
   more eigenvalues=  21.1483  34.5296  53.2843  
	 iter =   56     E =     1.80366783    dE =     0.00401382 
	 iter =   57     E =     1.79642138    dE =     0.00403383 
	 iter =   58     E =     1.78920518    dE =     0.00403319 
	 iter =   59     E =     1.78365836    dE =     0.00310980 
	 iter =   60     E =     1.58733421    dE =     0.12368167 
   more eigenvalues=  21.0927  34.4563  53.2355  
	 iter =   61     E =     1.52907467    dE =     0.03810118 
	 iter =   62     E =     1.52499002    dE =     0.00267848 
	 iter =   63     E =     1.51444560    dE =     0.00696256 
	 iter =   64     E =     1.51042556    dE =     0.00266152 
	 iter =   65     E =     1.50733101    dE =     0.00205300 
   more eigenvalues=  20.7449  33.3795  45.4292  
	 iter =   66     E =     1.50217469    dE =     0.00343257 
	 iter =   67     E =     1.49627680    dE =     0.00394171 
	 iter =   68     E =     1.49360334    dE =     0.00178994 
	 iter =   69     E =     1.48528864    dE =     0.00559804 
	 iter =   70     E =     1.39563461    dE =     0.06423890 
   more eigenvalues=  20.7192  33.3244  45.424  
	 iter =   71     E =     1.16191313    dE =     0.20115228 
	 iter =   72     E =     1.16118954    dE =     0.00062315 
	 iter =   73     E =     1.15837200    dE =     0.00243232 
	 iter =   74     E =     1.15573266    dE =     0.00228370 
	 iter =   75     E =     1.15495642    dE =     0.00067209 
   more eigenvalues=  20.6896  33.2532  45.4014  
	 iter =   76     E =     1.15406940    dE =     0.00076860 
	 iter =   77     E =     1.13398709    dE =     0.01770947 
	 iter =   78     E =     1.13357987    dE =     0.00035923 
	 iter =   79     E =     1.13345271    dE =     0.00011219 
	 iter =   80     E =     1.13237723    dE =     0.00094975 
   more eigenvalues=  18.4398  30.0028  41.7954  
	 iter =   81     E =     1.13090239    dE =     0.00130413 
	 iter =   82     E =     1.13064943    dE =     0.00022373 
	 iter =   83     E =     1.12708797    dE =     0.00315988 
	 iter =   84     E =     1.12676881    dE =     0.00028325 
	 iter =   85     E =     1.12645672    dE =     0.00027706 
   more eigenvalues=  18.2448  28.8981  40.0402  
	 iter =   86     E =     1.12610721    dE =     0.00031037 
	 iter =   87     E =     1.12580421    dE =     0.00026914 
	 iter =   88     E =     1.12515339    dE =     0.00057843 
	 iter =   89     E =     1.11932130    dE =     0.00521038 
	 iter =   90     E =     1.11837920    dE =     0.00084238 
   more eigenvalues=  18.2385  28.8759  40.018  
	 iter =   91     E =     1.11829847    dE =     0.00007218 
	 iter =   92     E =     1.11818895    dE =     0.00009794 
	 iter =   93     E =     1.11795966    dE =     0.00020510 
	 iter =   94     E =     1.11773808    dE =     0.00019825 
	 iter =   95     E =     1.11764276    dE =     0.00008528 
   more eigenvalues=  18.2059  28.8161  39.9782  
	 iter =   96     E =     1.11734012    dE =     0.00027086 
	 iter =   97     E =     1.11724716    dE =     0.00008321 
	 iter =   98     E =     1.11696444    dE =     0.00025311 
	 iter =   99     E =     1.11631386    dE =     0.00058280 
	 iter =  100     E =     1.11619976    dE =     0.00010222 
   more eigenvalues=  18.205  28.8127  39.974  
	 iter =  101     E =     1.11604139    dE =     0.00014191 
	 iter =  102     E =     1.11563034    dE =     0.00036844 
	 iter =  103     E =     1.11546866    dE =     0.00014495 
	 iter =  104     E =     1.11506852    dE =     0.00035885 
finding new state with lower energy failed

time=  41.03
done!
пн фев 18 22:56:16 IST 2019
rm: cannot remove ‘sbatch_script.deuteron_Cut6_eB_0.7.16812’: No such file or directory
