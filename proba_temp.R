death_female < function(age, ethrisk, b_type1, b_type2, bmi,surv){
  survivor = c(rep(0,90), 0.999977290630341)

  #Ichemocat = c(0, 0.8345234517964427167768804, 1.2595202625466794810193960, 2.8514552483447963560081462)
  Iethrisk = c(0,0,
    0.6389369130057840351355480,
    0.3342933389380337572127644,
    0.3428258158976604796919219,
    0.1716307703741346002423995,
    0.5199930351630326352818656,
    0.6823609168626041387994974,
    0.1930811621976745162676536,
    0.5483388756920363205082936
  );
  dage = age/10;
  age_2 = dage^3*log(dage);
  age_1 = dage^3;
  dbmi = dbmi/10;
  bmi_1 = (dbmi)^(0.5);
  bmi_2 = (dbmi)^(0.5)*log(dbmi);
  age_1 = age_1 - 115.599884033203125;
  age_2 = age_2 - 183.038345336914062;
  bmi_1 = bmi_1 - 1.632479429244995;
  bmi_2 = bmi_2 - 1.600156426429749;
  #town = town - 0.327639788389206;
  
  a=0;
  a = a +  Iethrisk[ethrisk];
  a = a +  age_1 * 0.0535266800950749549459218;
  a = a +  age_2 * -0.0200935878258154260178614;
  a = a +  bmi_1 * -19.7435582245984164728724863;
  a = a +  bmi_2 * 6.6648702078668167203545636;
  #a = a +  town * 0.0787269477751315061020421;
  
  #a += b2_82 * 0.0859851843797995313289917;
  a = a + b_type1 * 1.3918300744178950800744587;
  a = a +  b_type2 * 1.8389652399973426266654997;
  a = a +  age_1 * b_type2 * -0.0200621605517602719093162;
  a = a +  age_2 * b_type2 * 0.0074957790032429043661222;
  
  score = 100.0 * (1 - survivor[surv]^(exp(a)) );
  return score;
}


/* 
  * QCOVID © Copyright, Oxford University 2020.
* All Rights Reserved. The author, being Professor Julia Hippsley-Cox, has asserted their moral right.
* 
  */
  
  static double hospital_female(
    int age,int b2_82,int b2_leukolaba,int b2_prednisolone,int b_AF,int b_CCF,int b_asthma,int b_bloodcancer,int b_cerebralpalsy,int b_chd,int b_cirrhosis,int b_congenheart,int b_copd,int b_dementia,int b_epilepsy,int b_fracture4,int b_neurorare,int b_parkinsons,int b_pulmhyper,int b_pulmrare,int b_pvd,int b_ra_sle,int b_respcancer,int b_semi,int b_sicklecelldisease,int b_stroke,int b_type1,int b_type2,int b_vte,double bmi,int chemocat,int ethrisk,int homecat,int learncat,int p_marrow6,int p_radio6,int p_solidtransplant,int renalcat,int surv,double town
  )
{
  double survivor[91] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.999614596366882
  };
  
  double Ichemocat[4] = {
    0,
    0.7468184809678632962715028,
    1.4335568622103449509808115,
    2.7425623418828077859643599
  };
  double Iethrisk[10] = {
    0,
    0,
    0.6381416181144948795989080,
    0.4154336125176559812999244,
    0.3460373664440339336323404,
    0.7629835412388623616664063,
    0.6956993818885366387405611,
    0.8317012026776886557399848,
    0.1388269618785758774404115,
    0.6432491845097410010367867
  };
  double Ihomecat[3] = {
    0,
    0.6111109399113108242573844,
    0.2046773133743576833509792
  };
  double Ilearncat[3] = {
    0,
    0.4267124890741003651051244,
    2.1792421003277180346913156
  };
  double Irenalcat[7] = {
    0,
    0,
    0.2982583643377000881535821,
    0.5848878422963580403504125,
    1.4268609433854653190110184,
    1.3148983941801617447708850,
    1.7128356862419007455855535
  };
  
  double dage = age;
  dage=dage/10;
  double age_1 = pow(dage,.5);
  double age_2 = pow(dage,2);
  double dbmi = bmi;
  dbmi=dbmi/10;
  double bmi_2 = log(dbmi);
  double bmi_1 = pow(dbmi,-2);
  
  age_1 = age_1 - 2.207121372222900;
  age_2 = age_2 - 23.730392456054688;
  bmi_1 = bmi_1 - 0.140802085399628;
  bmi_2 = bmi_2 - 0.980199992656708;
  town = town - 0.327639788389206;
  
  double a=0;
  
  a += Ichemocat[chemocat];
  a += Iethrisk[ethrisk];
  a += Ihomecat[homecat];
  a += Ilearncat[learncat];
  a += Irenalcat[renalcat];
  
  a += age_1 * -0.1484733673762321515265938;
  a += age_2 * 0.0405941535676193412940371;
  a += bmi_1 * 6.1144623880326562925802136;
  a += bmi_2 * 2.7351660262730592698687815;
  a += town * 0.0837552324383479818159515;
  
  a += b2_82 * 0.2742900245435872519372822;
  a += b2_leukolaba * 0.2676801383148893487273767;
  a += b2_prednisolone * 0.6545953126312467063030454;
  a += b_AF * 0.2903016741277205103877179;
  a += b_CCF * 0.3236173092786057137182354;
  a += b_asthma * 0.1121787272583002481596282;
  a += b_bloodcancer * 0.3377863143614412422266469;
  a += b_cerebralpalsy * 0.9778476962006689143791505;
  a += b_chd * 0.1076997815345036024758940;
  a += b_cirrhosis * 0.6047183940676497115873644;
  a += b_congenheart * 0.1825617255143640038639319;
  a += b_copd * 0.2940790450123081933853086;
  a += b_dementia * 0.5469427350294541190223185;
  a += b_epilepsy * 0.4527175076655778340750658;
  a += b_fracture4 * 0.2983585791165178080497355;
  a += b_neurorare * 0.9059625435491620581984762;
  a += b_parkinsons * 0.5286796359462073713331165;
  a += b_pulmhyper * 0.4692582807440949799193675;
  a += b_pulmrare * 0.2499579320244159075237178;
  a += b_pvd * 0.1942006598368692937839342;
  a += b_ra_sle * 0.3019421121979438682458863;
  a += b_respcancer * 0.4999823548811314077866541;
  a += b_semi * 0.3162439787347626762858965;
  a += b_sicklecelldisease * 1.8985404911409835548852243;
  a += b_stroke * 0.3305273237679022813040319;
  a += b_type1 * 1.3943531760373193417734683;
  a += b_type2 * 0.9708326823437052333076736;
  a += b_vte * 0.2955542808178847624667185;
  a += p_marrow6 * 0.1893906144867409657717161;
  a += p_radio6 * 0.3876198174592248579806153;
  a += p_solidtransplant * 0.4490079929764763111421644;
  
  a += age_1 * b_type2 * -1.1514860942738034399468461;
  a += age_2 * b_type2 * 0.0018396028070442396740169;
  
  double score = 100.0 * (1 - pow(survivor[surv], exp(a)) );
  return score;
  }