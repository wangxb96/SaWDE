function [gbestPos,gbestObjVal,outGeEvaNumGbest,outStrSelPRecord]= LSFS(fhd,Max_Gen,Max_FES,Particle_Number,Dimension,vecNum,dataset,pro)
% We are more than happy to share the related Matlab codes with you. If you do some further studies based on this research work,  please cite this study as a reference in your paper(s) as follows:
%  Yu Xue, Bing Xue, and Mengjie Zhang. 2019. Self-adaptive Particle Swarm Optimization for Large-scale Feature Selection in Classification. ACM Transactions on Knowledge Discovery from Data. 1, 1, Article 1 (January 2019), 26 pages. https://doi.org/10.1145/3340848
% Please notice that the reference format and the codes are not the final versions, we will continuously update them in the following some days. Please feel free to contact me if you have any questions or comments. My Email address is: xueyu@nuist.edu.cn
% I will also update the reference and codes in my ResearchGate website at 
% https://www.researchgate.net/profile/Yu_Xue19?ev=hdr_xprf
% 
% (1) File "testingFrameworkOnXue" is the main procedure
% (2) The maximum number of fitness evaluation can be modified in the file "\LSFS\algLib\Benchmarks\Xue\getFunParamOnXue.m"
% (3) The codes for the proposed method are included in the file of "\LSFS\algLib\algorithms\LSFS"
%Setting parameter values

threshold=0.6;

GeEvaNumGbest_SaPSOLSFS=[];    


[XRmin1, XRmax1, Lbound1, Ubound1, oo, maxfeval1] = getFunParamOnXue( Dimension, pro );


%-------------------------------------------
% step1 Initialization 
%-------------------------------------------
rand('state',sum(100*clock));  
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
funSeq =pro;  
fitcount=0;


%##########################################################################

 [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS] = initialPopulation01(ps,D);  
 [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,fitcount,...   
    pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,gbestrep_SaPSOLSFS...  
    ] = ...    
    evalPop_SaPSOLSFS(fhd,dataset,funSeq,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,ps,D,fitcount);  


[val_SaPSOLSFS,ind_SaPSOLSFS] = sort(popVal_SaPSOLSFS);  


straNum_SaPSOLSFS = 5;
boundStraNum_SaPSOLSFS=straNum_SaPSOLSFS+2;  
SSP_SaPSOLSFS=zeros(1,straNum_SaPSOLSFS);     
SSP_SaPSOLSFS(1,:)=1/straNum_SaPSOLSFS;       
vSSP_SaPSOLSFS=repmat(SSP_SaPSOLSFS,ps,1);   
 
LP_SaPSOLSFS = randSelNum(10);             
 

 nsFlagSaPSOLSFS=zeros(ps,straNum_SaPSOLSFS);           
 nfFlagSaPSOLSFS=zeros(ps,straNum_SaPSOLSFS);
 iniSSP_SaPSOLSFS = SSP_SaPSOLSFS;                  
 

 CRm=zeros(1,straNum_SaPSOLSFS);     
 CRm(1,:)=0.5;
 CRMemory=zeros(ps,straNum_SaPSOLSFS); 
 
 Skg=[];       
 Fkg=[];
 TempCRMem=[]; 
 

eachTopNum_SaPSOLSFS=3;

ringTop_SaPSOLSFS = genRingToplogy(ps,eachTopNum_SaPSOLSFS);         

curIter=0;

optimal_solution=find(gbest_SaPSOLSFS>=threshold);
optimal_solution_size=size(optimal_solution,2);


outCell={curIter,fitcount-fitcount,optimal_solution,optimal_solution_size,gbestval_SaPSOLSFS};
GeEvaNumGbest_SaPSOLSFS=[GeEvaNumGbest_SaPSOLSFS;outCell];


%-------------------------------------------%
% step2 
%-------------------------------------------%
                              
fc_SaPSOLSFS=[];                   
straSelPro_SaPSOLSFS = [];         
curIterLast=1;                
curIter=1;                    
lastIterLP_SaPSOLSFS=0;          


      

while fitcount<=Max_FES %main loop
    
    curIter = curIter+1;
%-------------------------------------------------------------------------

    
    improveInd_SaPSOLSFS=[];           
    for curParticle=1:ps            
               
            strategySeq = getStrategyOnAParticle(vSSP_SaPSOLSFS,curParticle,ringTop_SaPSOLSFS,3*straNum_SaPSOLSFS);  
                        
            switch strategySeq    

                  case 1
                   [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS16(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS);
                  case 2
                   [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS18(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS);
                  case 3
                   [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS17(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS);
                  case 4                
                   [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS2(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS);
                  case 5
                   [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS14(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS);

                otherwise    
                  %do somethings 
            end
%                   straPrint=sprintf('strategy %02d is selected ',strategySeq);
%                   disp(straPrint);
            
if mod(fitcount,10000)==0 %store the information outputs
fitcount
optimal_solution=find(gbest_SaPSOLSFS>=threshold);
optimal_solution_size=size(optimal_solution,2)
gbestval_SaPSOLSFS
outCell={curIter,fitcount,optimal_solution,optimal_solution_size,gbestval_SaPSOLSFS};

GeEvaNumGbest_SaPSOLSFS=[GeEvaNumGbest_SaPSOLSFS;outCell];
end
             SSPFlag_SaPSOLSFS(curParticle,strategySeq) = 1;


end    
    

    
    [val_SaPSOLSFS,ind_SaPSOLSFS]=sort(popVal_SaPSOLSFS);
   
%     [TempCRMem,CRMemory]=uppdateCRMeInOneGen(nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,ps,straNum_SaPSOLSFS,TempCRMem); 
    [Skg,Fkg,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS]=uppdateNsNfInOneGen(nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,straNum_SaPSOLSFS,Skg,Fkg,ps); 

 
    if curIter-lastIterLP_SaPSOLSFS==LP_SaPSOLSFS
 
        lastIterLP_SaPSOLSFS = curIter;
        [SSP_SaPSOLSFS,Skg,Fkg] = updateSSP_SaPSOLSFS(Skg,Fkg,SSP_SaPSOLSFS,straNum_SaPSOLSFS);    
        outSSP=SSP_SaPSOLSFS;
        straSelPro_SaPSOLSFS = [straSelPro_SaPSOLSFS; outSSP];
        fc_SaPSOLSFS = [fc_SaPSOLSFS fitcount];
        vSSP_SaPSOLSFS=repmat(SSP_SaPSOLSFS,ps,1);
        
%         [CRm,TempCRMem]=uppdateCRm(CRm,TempCRMem); 
    end
    


end 


StrSelPRecord_SaPSOLSFS=zeros(size(straSelPro_SaPSOLSFS,1),boundStraNum_SaPSOLSFS);
for i=1:size(straSelPro_SaPSOLSFS,1) 

     StrSelPRecord_SaPSOLSFS(i,1)=i;           
     StrSelPRecord_SaPSOLSFS(i,2)=fc_SaPSOLSFS(i); 
     StrSelPRecord_SaPSOLSFS(i,3:boundStraNum_SaPSOLSFS)=straSelPro_SaPSOLSFS(i,:);  
end



gbestPos= optimal_solution;
gbestObjVal=gbestval_SaPSOLSFS;



    outGeEvaNumGbest=GeEvaNumGbest_SaPSOLSFS;

    outStrSelPRecord=StrSelPRecord_SaPSOLSFS;

  

function ringTop = genRingToplogy(ps,eachTopNum)

ringTop = zeros(ps,eachTopNum);
for i=1:ps
    if i>ps-eachTopNum+1
        ringTop(i,1:(ps-i+1)) = i:ps;        
        ringTop(i,(ps-i+1+1):eachTopNum) = 1:(eachTopNum-(ps-i+1));
    else
        ringTop(i,:) = i:(i+eachTopNum-1);
    end
end


function strategySeq = getStrategyOnAParticle(SSP_CSA,curParticle,ringTop,num)


curParticleTop = ringTop(curParticle,:);    
n = size(curParticleTop,2);                 
P = sum(SSP_CSA(curParticleTop(1:n),:),1)./n;  
strategySeq = roulette(P,num);            

function strategySeq = roulette(P,num)


m = length(P); 
Select = zeros(1,num); 
r = rand(1,num); 
for i=1:num
    sumP = 0;
    r2 = randperm(m);
    j = r2(1); 
    while sumP < r(i)
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end

    Select(i) = mod(j-2,m)+1;
end

tempSelect = zeros(1,m);
for i=1:m
    tempSelect(i) = sum(Select==i);
end
[~,ind_CSA] = sort(tempSelect,'descend');
strategySeq = ind_CSA(1);



function [VRminMat,VRmaxMat] = initialCommParam(ps,D,VRmin,VRmax)


if length(VRmin)==1            
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.2*(VRmax-VRmin);           
Vmin=repmat(-mv,ps,1);          
Vmax=-Vmin;
VRminMat=repmat(VRmin,ps,1);    
VRmaxMat=repmat(VRmax,ps,1);

function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS]= initialPopulation01(ps,D)

population_SaPSOLSFS_conti=rand(ps,D);       
population_SaPSOLSFS_veloci=rand(ps,D)-0.5; 
p1=0.8;
p2=0.8;
population_SaPSOLSFS=zeros(ps,D);  
for i=1:ps
    r1=rand;
    if r1<=p1  
        for j=1:D  
            if population_SaPSOLSFS_conti(i,j)>=p2
                population_SaPSOLSFS(i,j)=1;
             else
                population_SaPSOLSFS(i,j)=0;
            end
        end
        if size(find(population_SaPSOLSFS(i,:)==0),2)==D  
        r3=randperm(D);
        r4=r3(1);
        population_SaPSOLSFS(i,r4)=1;
        population_SaPSOLSFS_conti(i,r4)=p2+(1-p2)*rand;
        end
        if size(find(population_SaPSOLSFS(i,:)==1),2)==D   
        r3=randperm(D);
        r4=r3(1);
        population_SaPSOLSFS(i,r4)=0;
        population_SaPSOLSFS_conti(i,r4)=p2*rand;
        end
    else      
        for j=1:D  
            if population_SaPSOLSFS_conti(i,j)<=p2
                population_SaPSOLSFS(i,j)=1;
                population_SaPSOLSFS_conti(i,j)=p2+(1-p2)*rand;
            else
                population_SaPSOLSFS(i,j)=0;
                population_SaPSOLSFS_conti(i,j)=p2*rand;
            end
        end
        if size(find(population_SaPSOLSFS(i,:)==0),2)==D 
        r3=randperm(D);
        r4=r3(1);
        population_SaPSOLSFS(i,r4)=1;
        population_SaPSOLSFS_conti(i,r4)=p2+(1-p2)*rand;
        end
        if size(find(population_SaPSOLSFS(i,:)==1),2)==D  
        r3=randperm(D);
        r4=r3(1);
        population_SaPSOLSFS(i,r4)=0;
        population_SaPSOLSFS_conti(i,r4)=p2*rand;
        end
        
    end
end




function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,populationVal,fitcount,... 
    pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,gbestrep_SaPSOLSFS... 
    ] = ... %
    evalPop_SaPSOLSFS(fhd,dataset,funSeq,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,ps,D,fitcount)


for i=1:ps;    
    e(i,1)=feval(fhd,population_SaPSOLSFS(i,:),dataset,funSeq);  
end
fitcount=fitcount+ps;

populationVal =e;

pbest_SaPSOLSFS=population_SaPSOLSFS_conti;            
pbestval_SaPSOLSFS=e;                             
[gbestval_SaPSOLSFS,gbestid]=min(pbestval_SaPSOLSFS);  
gbest_SaPSOLSFS=population_SaPSOLSFS_conti(gbestid,:); 
gbestrep_SaPSOLSFS=repmat(gbest_SaPSOLSFS,ps,1);      


function [Vmin,Vmax,mv] = setVelLim(ps,D,VRmin,VRmax)


if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);        
Vmin=repmat(-mv,ps,1);      
Vmax=-Vmin;

function Vold=initialVold(ps,D,Vmin,Vmax)
    
    Vold=Vmin+rand(ps,D).*(Vmax-Vmin);

function memSize = setMemParam(ps)


memSize = ceil(ps);

     


function r= randACauchy(varargin)


medianVal = 0;
upperQuartileVal = 1;
n = rand();

if(nargin >= 1)
    medianVal=	varargin{1};
    if(nargin >= 2)
        upperQuartileVal=			varargin{2};
        upperQuartileVal(upperQuartileVal <= 0)=	NaN;	
        if(nargin >= 3),	n=	[varargin{3:end}];		
        end
    end
end
r = upperQuartileVal^2/(pi*(n^2+upperQuartileVal^2));



function new = ctrlVelLim(old,mv)


new=(old>mv).*mv+(old<=mv).*old; 

function newAg = ctrlVarBoundary(oldAg,VRmin,VRmax)

v1 = oldAg;
lu(1,:) = VRmin;
lu(2,:) = VRmax;
vioLow = find(v1 < lu(1, :));
if ~isempty(vioLow)
    v1(1, vioLow) = 2 .* lu(1, vioLow) - v1(1, vioLow);
    vioLowUpper = find(v1(1, vioLow) > lu(2, vioLow));
    if ~isempty(vioLowUpper)
        v1(1, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
    end
end
vioUpper = find(v1 > lu(2, :));
if ~isempty(vioUpper)
    v1(1, vioUpper) = 2 .* lu(2, vioUpper) - v1(1, vioUpper);
    vioUpperLow = find(v1(1, vioUpper) < lu(1, vioUpper));
    if ~isempty(vioUpperLow)
        v1(1, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
    end
end
newAg=v1;


function SRR = getSecondResRate()

SRR = randSelNum([0.1,0.3,0.9]);


function re = randANum(lo,up)

re = randNum(lo,up,1,1);

function re = randNum(lo,up,m,n)

re =lo+(up-lo).*rand(m,n);

function re = randSelNum(NumList)

r = randperm(length(NumList));
re = NumList(r(1));




function [Skg,Fkg,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS]=uppdateNsNfInOneGen(nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,straNum_SaPSOLSFS,Skg,Fkg,ps)
TempSkg=sum(nsFlagSaPSOLSFS,1);
Skg=[Skg;TempSkg];
TempFkg=sum(nfFlagSaPSOLSFS,1);
Fkg=[Fkg;TempFkg];
nsFlagSaPSOLSFS=zeros(ps,straNum_SaPSOLSFS);
nfFlagSaPSOLSFS=zeros(ps,straNum_SaPSOLSFS);

function [SSP_SaPSOLSFS,Skg,Fkg]=updateSSP_SaPSOLSFS(Skg,Fkg,SSP_SaPSOLSFS,straNum_SaPSOLSFS)
        tempSum1=sum(Skg,1);    
        tempSum2=tempSum1;
        [i,j]=size(tempSum1);
        for k=1:j
           if (tempSum1(k)==0)
              tempSum2(k)=1;
           end
        end
        if (sum(tempSum1(1,:),2)==0)
        SSP_SaPSOLSFS(1,:)=1/straNum_SaPSOLSFS;
        else
        SSP_SaPSOLSFS=tempSum1./(tempSum2+sum(Fkg,1));
        SSP_SaPSOLSFS=SSP_SaPSOLSFS./sum(SSP_SaPSOLSFS,2);
        end
        Skg=[];
        Fkg=[];
        

        
        

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% The used CSGSs
function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS2(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS)
w=0.7298; 
c1=1.49618;
c2=1.49618;
r1=rand;
r2=rand;
threshold=0.6;
v_new=w.*population_SaPSOLSFS_veloci(curParticle,:)+c1*r1.*(pbest_SaPSOLSFS(curParticle,:)-population_SaPSOLSFS_conti(curParticle,:))+c2*r2.*(gbest_SaPSOLSFS-population_SaPSOLSFS_conti(curParticle,:)); 
population_SaPSOLSFS_veloci(curParticle,:)=v_new;
x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+v_new;  

ind1=find(x_new_conti<=0); 
if ~isempty(ind1)
x_new_conti(ind1)=rand(1,size(ind1,2));
end
ind1=find(x_new_conti>=1); 
if ~isempty(ind1)
x_new_conti(ind1)=rand(1,size(ind1,2));
end

num=0;
while size(find(x_new_conti(1,:)>=threshold),2)==0 
    num=num+1;
    if mod(num,50)==0
    population_SaPSOLSFS_conti(curParticle,:)=rand(1,D);
    end
    if num==2000
    sprintf('always sycle')
    pause;  
    end
    v_new=w.*population_SaPSOLSFS_veloci(curParticle,:)+c1*r1.*(pbest_SaPSOLSFS(curParticle,:)-population_SaPSOLSFS_conti(curParticle,:))+c2*r2.*(gbest_SaPSOLSFS-population_SaPSOLSFS_conti(curParticle,:));
    population_SaPSOLSFS_veloci(curParticle,:)=v_new;
    x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+v_new; 
   
   ind1=find(x_new_conti<=0);
   if ~isempty(ind1)
   x_new_conti(ind1)=rand(1,size(ind1,2));
   end
   ind1=find(x_new_conti>=1); 
   if ~isempty(ind1)
   x_new_conti(ind1)=rand(1,size(ind1,2));
   end

end

x_new_binary=zeros(1,D);  
ind=find(x_new_conti(1,:)>=threshold); 
x_new_binary(1,ind)=1;

x_new_binary_value=feval(fhd,x_new_binary,dataset,funSeq);  
fitcount=fitcount+1;


if (x_new_binary_value<popVal_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==popVal_SaPSOLSFS(curParticle)) && (size(ind,2)<size(find(population_SaPSOLSFS_conti(curParticle,:)>=threshold),2)))...
        ||(x_new_binary_value<pbestval_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2)))
   nsFlagSaPSOLSFS(curParticle,strategySeq)=1;
   improveInd_SaPSOLSFS=[improveInd_SaPSOLSFS,curParticle];
else
   nfFlagSaPSOLSFS(curParticle,strategySeq)=1;
end  



if x_new_binary_value<pbestval_SaPSOLSFS(curParticle)   
    pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
    pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;
end

if (x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2))
   pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
   pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;    
end
   

if x_new_binary_value<gbestval_SaPSOLSFS 
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;
end
if (x_new_binary_value==gbestval_SaPSOLSFS )&&(size(ind,2)<size(find(gbest_SaPSOLSFS>=threshold),2))
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;   
end

population_SaPSOLSFS_conti(curParticle,:)=x_new_conti; 
population_SaPSOLSFS(curParticle,:)=x_new_binary;          
popVal_SaPSOLSFS(curParticle)=x_new_binary_value;      



function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS14(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS)
threshold=0.6;

Vmin=-0.5.*ones(1,D);
Vmax=0.5.*ones(1,D);
Pmin=zeros(1,D);
Pmax=ones(1,D);

r1=rand(1,D);
r2=rand(1,D);
r3=rand(1,D);
r4=randperm(ps);

x_new_conti=r1.*population_SaPSOLSFS_conti(curParticle,:)+r2.*gbest_SaPSOLSFS+r3.*(population_SaPSOLSFS_conti(r4(1),:)-population_SaPSOLSFS_conti(r4(2),:)); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);

num=0;
while size(find(x_new_conti(1,:)>=threshold),2)==0 
    num=num+1;
    if mod(num,50)==0
    population_SaPSOLSFS_conti(curParticle,:)=rand(1,D);
    end
    if num==2000
    sprintf('always sycle')
    pause;  
    end

r1=rand(1,D);
r2=rand(1,D);
r3=rand(1,D);
r4=randperm(ps);
x_new_conti=r1.*population_SaPSOLSFS_conti(curParticle,:)+r2.*pbest_SaPSOLSFS(curParticle,:)+r3.*(population_SaPSOLSFS_conti(r4(1),:)-population_SaPSOLSFS_conti(r4(2),:));   
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);    


end

x_new_binary=zeros(1,D);  
ind=find(x_new_conti(1,:)>=threshold); 
x_new_binary(1,ind)=1;

x_new_binary_value=feval(fhd,x_new_binary,dataset,funSeq); 
fitcount=fitcount+1;


if (x_new_binary_value<popVal_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==popVal_SaPSOLSFS(curParticle)) && (size(ind,2)<size(find(population_SaPSOLSFS_conti(curParticle,:)>=threshold),2)))...
        ||(x_new_binary_value<pbestval_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2)))
   nsFlagSaPSOLSFS(curParticle,strategySeq)=1;
   improveInd_SaPSOLSFS=[improveInd_SaPSOLSFS,curParticle];
else
   nfFlagSaPSOLSFS(curParticle,strategySeq)=1;
end  



if x_new_binary_value<pbestval_SaPSOLSFS(curParticle)   
    pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
    pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;
end

if (x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2))
   pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
   pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;    
end
   

if x_new_binary_value<gbestval_SaPSOLSFS  
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;
end
if (x_new_binary_value==gbestval_SaPSOLSFS )&&(size(ind,2)<size(find(gbest_SaPSOLSFS>=threshold),2))
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;   
end

population_SaPSOLSFS_conti(curParticle,:)=x_new_conti; 
population_SaPSOLSFS(curParticle,:)=x_new_binary;          
popVal_SaPSOLSFS(curParticle)=x_new_binary_value;      



function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS16(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS)

threshold=0.6;
ebsilon=0.000000001;
c=(D-1)/D*normrnd(0,1,1,D)+1/D*randACauchy(0,1);
r=randperm(ps);

Vmin=-0.5.*ones(1,D);
Vmax=0.5.*ones(1,D);
Pmin=zeros(1,D);
Pmax=ones(1,D);

[val_popVal_SaPSOLSFS,ind_popVal_SaPSOLSFS]=sort(popVal_SaPSOLSFS);
sizePopulation=size(ind_popVal_SaPSOLSFS,1);
subSizePopulation=ceil(sizePopulation*0.5);
subPop=population_SaPSOLSFS_conti(ind_popVal_SaPSOLSFS(1:subSizePopulation),:);

meanVal=mean(subPop,1);


population_SaPSOLSFS_veloci(curParticle,:)=(meanVal(1,:)-population_SaPSOLSFS_conti(curParticle,:))+c/sqrt(3).*sqrt((pbest_SaPSOLSFS(curParticle,:)-meanVal(1,:)).^2+(population_SaPSOLSFS_conti(curParticle,:)-meanVal(1,:)).^2+(population_SaPSOLSFS_conti(r(1),:)-meanVal(1,:)).^2);; %����������ӵĲ������
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax); 

x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax); 

num=0;
while size(find(x_new_conti(1,:)>=threshold),2)==0 
    num=num+1;
    if mod(num,50)==0
    population_SaPSOLSFS_conti(curParticle,:)=rand(1,D);
    end
    if num==2000
    sprintf('always sycle')
    pause; 
    end

c=(D-1)/D*normrnd(0,1,1,D)+1/D*randACauchy(0,1);  
r=randperm(ps);


population_SaPSOLSFS_veloci(curParticle,:)=(meanVal(1,:)-population_SaPSOLSFS_conti(curParticle,:))+c/sqrt(3).*sqrt((pbest_SaPSOLSFS(curParticle,:)-meanVal(1,:)).^2+(population_SaPSOLSFS_conti(curParticle,:)-meanVal(1,:)).^2+(population_SaPSOLSFS_conti(r(1),:)-meanVal(1,:)).^2);; %����������ӵĲ������
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax); 


x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax); 
  
end

x_new_binary=zeros(1,D);  
ind=find(x_new_conti(1,:)>=threshold); 
x_new_binary(1,ind)=1;

x_new_binary_value=feval(fhd,x_new_binary,dataset,funSeq); 
fitcount=fitcount+1;


if (x_new_binary_value<popVal_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==popVal_SaPSOLSFS(curParticle)) && (size(ind,2)<size(find(population_SaPSOLSFS_conti(curParticle,:)>=threshold),2)))...
        ||(x_new_binary_value<pbestval_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2)))
   nsFlagSaPSOLSFS(curParticle,strategySeq)=1;
   improveInd_SaPSOLSFS=[improveInd_SaPSOLSFS,curParticle];
else
   nfFlagSaPSOLSFS(curParticle,strategySeq)=1;
end  



if x_new_binary_value<pbestval_SaPSOLSFS(curParticle)   
    pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
    pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;
end

if (x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2))
   pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
   pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;    
end
   

if x_new_binary_value<gbestval_SaPSOLSFS 
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;
end
if (x_new_binary_value==gbestval_SaPSOLSFS )&&(size(ind,2)<size(find(gbest_SaPSOLSFS>=threshold),2))
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;   
end

population_SaPSOLSFS_conti(curParticle,:)=x_new_conti;
population_SaPSOLSFS(curParticle,:)=x_new_binary;          
popVal_SaPSOLSFS(curParticle)=x_new_binary_value;      





function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS17(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS)


threshold=0.6;
ebsilon=0.000000001;
Vmin=-0.5.*ones(1,D);
Vmax=0.5.*ones(1,D);
Pmin=zeros(1,D);
Pmax=ones(1,D);

w=0.9-0.5*fitcount/Max_FES;
c=1.49445;
pbestf=zeros(1,D);
for i=1:D
    r=randperm(ps);
    pbestf(i)=pbest_SaPSOLSFS(r(1),i);
end

population_SaPSOLSFS_veloci(curParticle,:)=w.*population_SaPSOLSFS_veloci(curParticle,:)+c.*rand(1,D).*(pbestf(1,:)-population_SaPSOLSFS_conti(curParticle,:));
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax); 


x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);  

num=0;
while size(find(x_new_conti(1,:)>=threshold),2)==0 
    num=num+1;
    if mod(num,50)==0
    population_SaPSOLSFS_conti(curParticle,:)=rand(1,D);
    end
    if num==2000
    sprintf('always sycle')
    pause; 
    end

pbestf=zeros(1,D);
for i=1:D
    r=randperm(ps);
    pbestf(i)=pbest_SaPSOLSFS(r(1),i);
end

population_SaPSOLSFS_veloci(curParticle,:)=w.*population_SaPSOLSFS_veloci(curParticle,:)+c.*rand(1,D).*(pbestf(1,:)-population_SaPSOLSFS_conti(curParticle,:));
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax); 


x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);  
  
end

x_new_binary=zeros(1,D); 
ind=find(x_new_conti(1,:)>=threshold); 
x_new_binary(1,ind)=1;

x_new_binary_value=feval(fhd,x_new_binary,dataset,funSeq);  
fitcount=fitcount+1;


if (x_new_binary_value<popVal_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==popVal_SaPSOLSFS(curParticle)) && (size(ind,2)<size(find(population_SaPSOLSFS_conti(curParticle,:)>=threshold),2)))...
        ||(x_new_binary_value<pbestval_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2)))
   nsFlagSaPSOLSFS(curParticle,strategySeq)=1;
   improveInd_SaPSOLSFS=[improveInd_SaPSOLSFS,curParticle];
else
   nfFlagSaPSOLSFS(curParticle,strategySeq)=1;
end  



if x_new_binary_value<pbestval_SaPSOLSFS(curParticle)   
    pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
    pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;
end

if (x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2))
   pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
   pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;    
end
   

if x_new_binary_value<gbestval_SaPSOLSFS  
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;
end
if (x_new_binary_value==gbestval_SaPSOLSFS )&&(size(ind,2)<size(find(gbest_SaPSOLSFS>=threshold),2))
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;   
end

population_SaPSOLSFS_conti(curParticle,:)=x_new_conti; 
population_SaPSOLSFS(curParticle,:)=x_new_binary;          
popVal_SaPSOLSFS(curParticle)=x_new_binary_value;      




function [population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,fitcount,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS] = evolveWithStragegySaPSOLSFS18(fhd,dataset,funSeq,Max_Gen,Max_FES,fitcount,population_SaPSOLSFS_conti,population_SaPSOLSFS_veloci,population_SaPSOLSFS,popVal_SaPSOLSFS,pbest_SaPSOLSFS,pbestval_SaPSOLSFS,gbest_SaPSOLSFS,gbestval_SaPSOLSFS,ps,D,curParticle,CRm,nsFlagSaPSOLSFS,nfFlagSaPSOLSFS,CRMemory,improveInd_SaPSOLSFS,strategySeq,ringTop_SaPSOLSFS,eachTopNum_SaPSOLSFS)


threshold=0.6;
ebsilon=0.000000001;
Vmin=-0.5.*ones(1,D);
Vmax=0.5.*ones(1,D);
Pmin=zeros(1,D);
Pmax=ones(1,D);

w=0.9-0.5*fitcount/Max_FES;
c=1.49445;
pbestf=zeros(1,D);
for i=1:D
    r=randperm(ps);
    pbestf(i)=pbest_SaPSOLSFS(r(1),i);
end

population_SaPSOLSFS_veloci(curParticle,:)=w.*population_SaPSOLSFS_veloci(curParticle,:)+0.5*c.*rand(1,D).*(pbestf(1,:)-population_SaPSOLSFS_conti(curParticle,:)+pbest_SaPSOLSFS(curParticle,:)-population_SaPSOLSFS_conti(curParticle,:));
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax); 


x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);  

num=0;
while size(find(x_new_conti(1,:)>=threshold),2)==0 
    num=num+1;
    if mod(num,50)==0
    population_SaPSOLSFS_conti(curParticle,:)=rand(1,D);
    end
    if num==2000
    sprintf('always sycle')
    pause; 
    end

pbestf=zeros(1,D);
for i=1:D
    r=randperm(ps);
    pbestf(i)=pbest_SaPSOLSFS(r(1),i);
end

population_SaPSOLSFS_veloci(curParticle,:)=w.*population_SaPSOLSFS_veloci(curParticle,:)+0.5*c.*rand(1,D).*(pbestf(1,:)-population_SaPSOLSFS_conti(curParticle,:));
population_SaPSOLSFS_veloci(curParticle,:)=ctrlVarBoundary(population_SaPSOLSFS_veloci(curParticle,:),Vmin,Vmax);


x_new_conti=population_SaPSOLSFS_conti(curParticle,:)+population_SaPSOLSFS_veloci(curParticle,:); 
x_new_conti=ctrlVarBoundary(x_new_conti,Pmin,Pmax);  
  
end

x_new_binary=zeros(1,D);  
ind=find(x_new_conti(1,:)>=threshold); 
x_new_binary(1,ind)=1;

x_new_binary_value=feval(fhd,x_new_binary,dataset,funSeq); 
fitcount=fitcount+1;


if (x_new_binary_value<popVal_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==popVal_SaPSOLSFS(curParticle)) && (size(ind,2)<size(find(population_SaPSOLSFS_conti(curParticle,:)>=threshold),2)))...
        ||(x_new_binary_value<pbestval_SaPSOLSFS(curParticle))...
        ||((x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2)))
   nsFlagSaPSOLSFS(curParticle,strategySeq)=1;
   improveInd_SaPSOLSFS=[improveInd_SaPSOLSFS,curParticle];
else
   nfFlagSaPSOLSFS(curParticle,strategySeq)=1;
end  



if x_new_binary_value<pbestval_SaPSOLSFS(curParticle) 
    pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
    pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;
end

if (x_new_binary_value==pbestval_SaPSOLSFS(curParticle))&&(size(ind,2)<size(find(pbest_SaPSOLSFS(curParticle,:)>=threshold),2))
   pbest_SaPSOLSFS(curParticle,:)=x_new_conti;
   pbestval_SaPSOLSFS(curParticle)=x_new_binary_value;    
end
   

if x_new_binary_value<gbestval_SaPSOLSFS  
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;
end
if (x_new_binary_value==gbestval_SaPSOLSFS )&&(size(ind,2)<size(find(gbest_SaPSOLSFS>=threshold),2))
    gbest_SaPSOLSFS=x_new_conti;
    gbestval_SaPSOLSFS=x_new_binary_value;   
end

population_SaPSOLSFS_conti(curParticle,:)=x_new_conti; 
population_SaPSOLSFS(curParticle,:)=x_new_binary;          
popVal_SaPSOLSFS(curParticle)=x_new_binary_value;      




