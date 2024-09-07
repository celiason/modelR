#thickness of layers in unit cell
l1=0.2e-6
l2=0.8e-6

#permittivities of each layer
eps1=1
eps2=9

#lattice constant
a=l1+l2

#number of plane waves in Fourier series expansion
numG=10
N=numG*2+1

#counters for wave vector inside Brillouin zone
count_k=0
countG=1
countG1=1

chi<-data.frame(NULL)

for(G in seq(-numG*2*pi/a,numG*2*pi/a,2*pi/a))

{
	
	for(G1 in seq(-numG*2*pi/a,numG*2*pi/a,2*pi/a))

	{

	if((G-G1)==0)

	{ chi[countG1,countG]=1/(l1+l2)*(1/eps1*l1+1/eps2*l2) }
	
	else
	
{ 
	
chi[countG1,countG]=(0+1i)/(l1+l2)/(G-G1)*((1/eps1)*(exp((0-1i)*(G-G1)*l1)-1)+1/eps2*(exp((0-1i)*(G-G1)*(l1+l2))-exp((0-1i)*(G-G1)*l1))) 
	
}

countG=countG+1

}

countG1=countG1+1
countG=1

}

#reset loops
countG=1
countG1=1

#loop of band structure calculation

M<-data.frame(NULL)
dispe<-list(NULL)
ks<-seq(-pi/a,pi/a,.2*pi/a)

for(k in seq(-pi/a,pi/a,0.2*pi/a))

{
	for(G in seq(-numG*2*pi/a,numG*2*pi/a,2*pi/a))
	
	{
		
		for(G1 in seq(-numG*2*pi/a,numG*2*pi/a,2*pi/a))
		
		{
			
			M[countG1,countG]=chi[countG1,countG]*(k+G1)*(k+G)
			
			countG=countG+1
			
		}
			
			countG1=countG1+1
			countG=1
			
	}
			
			countG1=1
			
			V=eigen(M)$values
			count_k=count_k+1
			
			dispe[[count_k]]=sqrt((sort(abs(V))))*a/2/pi
			
}
plot(1,type='n',ylim=c(0,2),xlim=c(-pi/a,pi/a))
for(j in 1:8)

{

ys<-NULL

for(i in 1:11)

{

	temp<-dispe[[i]][j]
	ys<-c(ys,temp)
	
	}
lines(ys~ks,type='l',ylim=c(0,2),col='blue',lwd=2)}
	
