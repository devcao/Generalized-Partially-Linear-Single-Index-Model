A=read.table("Data Path")
dim(A)
B1=A[,8:27]
B1gt=A[,28:47]
B2=A[,48:67]
B_median=c()
B_median=c()
B_leftq=c()
B_rightq=c()
B_mean=c()
x_real=colSums(B2)/1000
B_real=x_real
B11=matrix(data=NA,nrow=1000,ncol=20)
for(i in 1:20)
{
 for(j in 1:1000)
{B11[j,i]=as.numeric(B1[j,i]+B1gt[j,i]*(x_real[i]-B2[j,i]))}
B_median[i]=median(B11[,i])
B_leftq[i]=quantile(B11[,i],0.025)
B_rightq[i]=quantile(B11[,i],0.975)
B_mean[i]=mean(B11[,i])
}
plot(x_real,B_real,type="l",col="black",ylim=c(-7.0,7.0),xlab='U',ylab='g')
lines(x_real,B_mean,col="red")
lines(x_real,B_median,col="green")
lines(x_real,B_leftq,col="blue")
lines(x_real,B_rightq,col="blue")

colSums(A[,2:7])/1000
median(A[,2])
sqrt(var(A[,2]))
median(A[,3])
sqrt(var(A[,3]))
median(A[,4])
sqrt(var(A[,4]))
median(A[,5])
sqrt(var(A[,5]))
median(A[,6])
sqrt(var(A[,6]))
median(A[,7])
sqrt(var(A[,7]))