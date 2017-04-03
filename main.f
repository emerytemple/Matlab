C	DESIGN OF INTERNAL-EXTERNAL EXPANSION PLUG NOZZLES
C	INPUT GAS=+1.0 WHEN DEALING WITH IDEAL GAS
C	INPUT GAS=-1.0 WHEN DEALING WITH REAL GAS
      DIMENSION HM(30),GAM(30)
	OPEN(20,file='input.dat')
	OPEN(30,file='output.dat')

      READ(20,*)N1,N2
C101   READ(5,*)N1,N2
C101   READ(5,11)N1,N2
11    FORMAT(2I4)

C      READ(5,1)GAS,PEIPC,XP,RRRE,RM,PHT
C      READ(5,*)GAS,PEIPC,XP,RRRE,RM,PHT
      READ(20,*)GAS,PEIPC,XP,RRRE,RM,PHT
C1     FORMAT(6F10.0)
C      WRITE (6,52) PEIPC
      WRITE (30,52) PEIPC
52    FORMAT(1H ,6HPEI/PC,1X,1H=,E14.7)
C      WRITE (6,53) XP
      WRITE (30,53) XP
53    FORMAT(1H ,9HEXPANSION,1X,5HRATIO,1X,1H=,F10.5)
C      WRITE (6,54) RRRE
      WRITE (30,54) RRRE
54    FORMAT(1H ,5HRR/RE,1X,1H=,E14.7)
C      WRITE (6,55) RM
      WRITE (30,55) RM
55    FORMAT(1H ,8HESTIMATE,1X,4HMACH,1X,6HNUMBER,1X,1H=,E14.7)
C      WRITE (6,56) PHT
      WRITE (30,56) PHT
56    FORMAT(1H ,3HPHT,1X,1H=,E14.7)

C      READ(5,66) R,TE,PAPC,G
C      READ(5,*) R,TE,PAPC,G
      READ(20,*) R,TE,PAPC,G
66    FORMAT(4F10.0)
C      WRITE(6,67) TE
      WRITE(30,67) TE
67    FORMAT(1H ,4HEXIT,1X,11HTEMPERATURE,1X,1H=,E14.7)
C      WRITE(6,68) PAPC
      WRITE(30,68) PAPC
68    FORMAT(1H ,5HPA/PC,1X,1H=,E14.7)
      IF(GAS)4,4,2

C2     READ(5,3) GAMA
C2     READ(5,*) GAMA
2     READ(20,*) GAMA
3     FORMAT(F5.0)
C      WRITE (6,57)
      WRITE (30,57)
57    FORMAT(1H ,5HUSING,1X,5HIDEAL,1X,3HGAS)
      GO TO 8

C4     READ(5,5) NT
C4     READ(5,*) NT
4     READ(20,*) NT
5     FORMAT(I4)
      DO 7 I=1,NT
C        READ(5,6) HM(I),GAM(I)
C        READ(5,*) HM(I),GAM(I)
        READ(20,*) HM(I),GAM(I)
6       FORMAT(2F10.7)
7     CONTINUE
C      WRITE(6,51)
      WRITE(30,51)
51    FORMAT(1H ,5HUSING,1X,4HREAL,1X,3HGAS)
8     CONTINUE
34    IF(GAS)30,9,9
30    DO 31 J=1,NT
        I=J
        IF(RM-HM(J))33,32,31
31    CONTINUE
32    GAMA=GAM(I)
      GO TO 9
33    GAMA=GAM(I-1)+(RM-HM(I-1))*(GAM(I)-GAM(I-1))/(HM(I)-HM(I-1))
9     FME=(2.0+(GAMA-1.0)*RM*RM)/(GAMA+1.0)
      COM=(GAMA+1.0)/(2.*(GAMA-1.0))
      FME=RM*XP-FME**COM
      FPM=(2.0+(GAMA-1.0)*RM*RM)/(GAMA+1.0)
      COM=(3.0-GAMA)/(2.0*(GAMA-1.0))
      FPM=XP-RM*(FPM**COM)
      DM=-FME/FPM
      RM=RM+DM
      DM=ABS(DM)
      IF(DM-0.00001) 10,10,34
10    CONTINUE 
      A=(GAMA-1.0)*(RM*RM-1.0)/(GAMA+1.0)
      A=SQRT(A)
      A=ATAN(A)
      B=SQRT((GAMA+1.0)/(GAMA-1.0))
      C=SQRT(RM*RM-1.0)
      C=ATAN(C)
      VE=B*A-C
      RMEI=PEIPC**((1.0-GAMA)/GAMA)
      RMEI=(2.0/(GAMA-1.))*(RMEI-1.0)
      RMEI=SQRT(RMEI)
      IF(GAS) 36,35,35
36    DO 37 J = 1, NT
        I=J
        IF(RMEI-HM(J))39,38,37
37    CONTINUE
38    GAMA=GAM(I)
      GO TO 35
39    GAMA=GAM(I-1)+(RM-HM(I-1))*(GAM(I)-GAM(I-1))/(HM(I)-HM(I-1))
35    A=(GAMA-1.0)*(RMEI*RMEI-1.0)/(GAMA+1.0)
      A=SQRT(A)
      A=ATAN(A)
      B=SQRT((GAMA+1.0)/(GAMA-1.0))
      C=SQRT(RMEI*RMEI-1.0)
      C=ATAN(C)
      VEI=B*A-C
      Y=1.0/RMEI
      UEI=ATAN(Y/SQRT(1.0-Y*Y))
      THEI=VEI-VE
      PHEI=THEI+UEI
C	CALCULATE THE ORIGIN OF THE LAST INTERNAL EXPANSION WAVE
      Al=(2.0/(GAMA+1.0))*(1.0+0.5*(GAMA-1.0)*RMEI*RMEI)
      A1=A1**((GAMA+1.0)/(2.0*(GAMA-1.0)))
      B1=SIN(PHEI)
      RPRE=SQRT(1.0-Al*B1/XP)
      XPRE=(RPRE-1.0)*COS(PHEI)/SIN(PHEI)
      XN1=N1
      DRM=(RMEI-1.0)/XN1
      K=0
      XM=1.0
C      WRITE(6,17)
      WRITE(30,17)
17    FORMAT(2H  ,4HMACH,11X,6HRX1/RE,9X,6HXX1/RE,9X,6HRX2/RE,9X,6HXX2/RE 
     1,9X,5HPX/PC)
44    IF(GAS)40,12,12
40    DO 43 J=1,NT
        I=J
        IF(XM-HM(J))41,42,43
43    CONTINUE
42    GAMA=GAM(I)
      GO TO 12
41    GAMA=GAM(I-1)+(XM-HM(I-1))*(GAM(I)-GAM(I-1))/(HM(I)-HM(I-1))
12    A=(GAMA-1.0)*(XM*XM-1.0)/(GAMA+1.0)
      A=SQRT(A) 
      A=ATAN(A) 
      B=SQRT((GAMA+1.0)/(GAMA-1.0))
      C=SQRT(XM*XM-1.0)
      C=ATAN(C)
      VX=A*B-C
      BX=PHT-1.570796-VX+ABS(THEI)
      XLRE=2.0*RRRE*SIN(0.5*BX)
      PSI=3.1416-PHT+VX-0.5*(3.1416-BX)
      RX1RE=RPRE+XLRE*SIN(PSI)
      XX1RE=XPRE-XLRE*COS(PSI)
      IF(K)62,60,62
60    RX2RE=SQRT(RX1RE*RX1RE+SIN(PHT)/XP)
      XX2RE=XX1RE+(RX2RE-RX1RE)*COS(PHT)/SIN(PHT)
      GO TO 61
62    UX=ATAN(1.0/(XM*SQRT(1.0-(1.0/XM)**2)))
      PHX=2.0*VEI-VE-VX+UX
      A2=(2.0/(GAMA+1.0))*(1.0+0.5*(GAMA-1.0)*XM*XM)
      B2=0.5*(GAMA+1.0)/(GAMA-1.0)
      RX2RE=SQRT(RX1RE*RX1RE+(A2**B2)*SIN(PHX)/XP)
      XX2RE=XX1RE+(RX2RE-RX1RE)*COS(PHX)/SIN(PHX)
61    PXPC=(1.0+0.5*(GAMA-1.0)*XM*XM)**(-GAMA/(GAMA-1.0))
C      WRITE(6,13) XM,RX1RE,XX1RE,RX2RE,XX2RE,PXPC
      WRITE(30,13) XM,RX1RE,XX1RE,RX2RE,XX2RE,PXPC
13    FORMAT(1H ,6(E14.7,1H ))
      K=K+1
      IF(K-N1) 14,14,15
14    XM=XM+DRM
      GO TO 44
15	WRITE(6,16)
  	WRITE(30,16)
16    FORMAT(1H ,8HINTERNAL,1X,7HPORTION,1X,2HOF,1X,3HTHE,1X,6HNOZZLE,1X
     1,2HIS,1X,8HCOMPLETE)
C	DESIGN OF EXTERNAL CONTOUR
C      WRITE(6,18)
      WRITE(30,18)
18    FORMAT(2H  ,4HMACH,10X,5HRX/RE,9X,5HXX/RE,9X,5HPX/PC,9X,5HCFVAC
     1,9X,3HSP.,1X,7HIMPULSE,3X,4HVAC.,1X,7HIMPULSE)
      UX=ATAN(1.0/(XM*SQRT(1.0-(1.0/XM)**2)))
      A=SQRT((GAMA-1.0)*(XM*XM-1.0)/(GAMA+1.0))
      A=ATAN(A)
      B=SQRT((GAMA+1.0)/(GAMA-1.0))
      C=SQRT(XM*XM-1.0)
      C=ATAN(C)
      VX=B*A-C
      RXRE=(2.0/(GAMA+1.0))*(1.0+0.5*(GAMA-1.0)*XM*XM)
      RXRE=RXRE**((GAMA+1.0)/(2.0*(GAMA-1.0)))
      RXRE=1.0-RXRE*SIN(VE-VX+UX)/XP
      RXRE=SQRT(RXRE)
      XXRE=(1.0-RXRE)*COS(VE-VX+UX)/SIN(VE-VX+UX)
      C1=(2.0/(GAMA+1.0))**(GAMA/(GAMA-1.0))
      C2=SQRT((0.5*(GAMA+1.0)*XM*XM)/(1.0+0.5*(GAMA-1.0)*XM*XM))
      SUMCG=GAMA*C1*C2*COS(THEI)+XP*PXPC*(1.0-RXRE*RXRE)
      VT=TE*(1.0+0.5*(GAMA-1.0)*RM*RM)
      VT=GAMA*G*R*VT/(0.5*(GAMA+1.0))
      VT=SQRT(VT)
      VQ=TE*(1.0+0.5*(GAMA-1.0)*RM*RM)
      VC=1.0+0.5*(GAMA-1.0)*XM*XM
      VQ=GAMA*R*G*VQ/VC
      VQ=SQRT(VQ)
      A=1.0+0.5*(GAMA-1.0)*XM*XM
      B=-GAMA/(GAMA-1.0)
      A=A**B
      C=1.0-PAPC*A*SIN(PHEI)/(XM*XM)
      D=COS(THEI)+C/GAMA
      SUMIM=VQ*D/G
      CO=COS(THEI)+1.0/GAMA
      SUMVA=VQ*CO/G
C      WRITE(6,19) XM,RXRE,XXRE,PXPC,SUMCG,SUMIM,SUMVA
      WRITE(30,19) XM,RXRE,XXRE,PXPC,SUMCG,SUMIM,SUMVA
19    FORMAT(1H ,7E14.7)
      K1=1
      XN2=N2
      DRM=(RM-XM)/XN2
      XM=XM+DRM
      PRO=PXPC
      RXO=RXRE
50    UX=ATAN (1.0/(XM*SQRT(1.0-(1.0/XM)**2)))
      IF(GAS)46,45,45
46    DO 49 J=1,NT
        I=J
        IF(XM-HM(J)) 47,48,49
49    CONTINUE
48    GAMA=GAM(I)
      GO TO 45
47    GAMA=GAM(I-1)+(XM-HM(I-1))*(GAM(9)-GAM(I-1))/(HM(I)-HM(I-1))
45    A=SQRT((GAMA-1.0)*(XM*XM-1.0)/(GAMA+1.0))
      A=ATAN(A)
      B=SQRT((GAMA+1.0)/(GAMA-1.0))
      C=SQRT(XM*XM-1.0)
      C=ATAN(C)
      VX=B*A-C
      RXRE=(2.0/(GAMA+1.0))*(1.0+0.5*(GAMA-1.0)*XM*XM)
      RXRE=RXRE**((GAMA+1.0)*0.5/(GAMA-1.0))
      RXRE=SQRT(1.0-RXRE*SIN(VE-VX+UX)/XP)
65    XXRE=(1.0-RXRE)*COS(VE-VX+UX)/SIN(VE-VX+UX)
      PXPC=(1.0+0.5*(GAMA-1.0)*XM*XM)**(-GAMA/(GAMA-1.0))
      SUMCG=SUMCG+0.5*XP*(PRO+PXPC)*(RXO*RXO-RXRE*RXRE)
      A=GAMA/(GAMA-1.0)
      A=(0.5*(GAMA+1.0))**A
      A=A*VT/(GAMA*G)
      B=0.5*A*XP
      SUMIM=SUMIM+B*(PRO+PXPC-2.0*PAPC)*(RXO*RXO-RXRE*RXRE)
      SUMVA=SUMVA+B*(PRO+PXPC)*(RXO*RXO-RXRE*RXRE)
C      WRITE(6,21) XM,RXRE,XXRE,PXPC,SUMCG,SUMIM,SUMVA
      WRITE(30,21) XM,RXRE,XXRE,PXPC,SUMCG,SUMIM,SUMVA
21    FORMAT(1H ,7E14.7)
      IF(K1-N2)22,23,23
22    XM=XM+DRM
      PRO=PXPC
      RXO=RXRE
      K1=K1+1
      GO TO 50
23    WRITE(6,24)
      WRITE(30,24)
24    FORMAT(1H ,8HEXTERNAL,1X,7HPORTION,1X,2HOF,1X,3HTHE,1X,6HNOZZLE,1X
     1,2HIS,1X,8HCOMPLETE)
C      GO TO 101
      END