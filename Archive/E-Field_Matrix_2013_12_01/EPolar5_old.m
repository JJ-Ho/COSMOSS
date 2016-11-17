function E5 = EPolar5(O_Sum,O_Vi,O_Probe,O_Pump2,O_Pump1)
%EPOLAR5
%    E5 = EPOLAR5(O_PUMP1,O_PUMP2,O_PROBE,O_SUM,O_VI)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    01-Dec-2013 08:13:59

t2 = cos(O_Pump2);
t3 = cos(O_Probe);
t4 = cos(O_Sum);
t5 = cos(O_Vi);
t6 = cos(O_Pump1);
t7 = sin(O_Pump1);
t8 = sin(O_Pump2);
t9 = sin(O_Probe);
t10 = sin(O_Vi);
t11 = sin(O_Sum);
E5 = [t2.*t3.*t4.*t5.*t6,t2.*t3.*t4.*t5.*t7,t3.*t4.*t5.*t6.*t8,t3.*t4.*t5.*t7.*t8,t2.*t4.*t5.*t6.*t9,t2.*t4.*t5.*t7.*t9,t4.*t5.*t6.*t8.*t9,t4.*t5.*t7.*t8.*t9,t2.*t3.*t4.*t6.*t10,t2.*t3.*t4.*t7.*t10,t3.*t4.*t6.*t8.*t10,t3.*t4.*t7.*t8.*t10,t2.*t4.*t6.*t9.*t10,t2.*t4.*t7.*t9.*t10,t4.*t6.*t8.*t9.*t10,t4.*t7.*t8.*t9.*t10,t2.*t3.*t5.*t6.*t11,t2.*t3.*t5.*t7.*t11,t3.*t5.*t6.*t8.*t11,t3.*t5.*t7.*t8.*t11,t2.*t5.*t6.*t9.*t11,t2.*t5.*t7.*t9.*t11,t5.*t6.*t8.*t9.*t11,t5.*t7.*t8.*t9.*t11,t2.*t3.*t6.*t10.*t11,t2.*t3.*t7.*t10.*t11,t3.*t6.*t8.*t10.*t11,t3.*t7.*t8.*t10.*t11,t2.*t6.*t9.*t10.*t11,t2.*t7.*t9.*t10.*t11,t6.*t8.*t9.*t10.*t11,t7.*t8.*t9.*t10.*t11];