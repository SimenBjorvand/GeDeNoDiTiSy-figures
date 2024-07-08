using PyPlot
using LinearAlgebra

#Function 
f1(p,x) = 4p[1]*x
f2(p,x) = 4*(1-p[1])*x + 2p[1] - 1
f3(p,x) = 4*(p[2]-1)*x + 3 - 2p[2]
f4(p,x) = 4*(1-x)*p[2]

f(p,x) = min(min(f1(p,x),f2(p,x)),max(f3(p,x),f4(p,x)))


# LD-derivative of function 
function SLmin(a::Array{Float64,1},b::Array{Float64,1})
    n = length(a)
    for i in 1:n
        if a[i] < b[i]
            return a[2:end]
        elseif b[i] < a[i]
            return b[2:end]
        end
    end
    return a[2:end]
end
function SLmin(A::Array{Float64,2},B::Array{Float64,2})
    r,c = size(A)
    M = Array{Float64,2}(undef,r,c-1)
    for i in 1:r
        M[i,:] = SLmin(A[i,:],B[i,:])
    end
    return M
end

Jf1(p,x) = [4x 0 4*p[1]]
Jf2(p,x) = [2-4x 0 4*(1-p[1])]
Jf3(p,x) = [0 4x-2 4*(p[2]-1)]
Jf4(p,x) = [0 4*(1-x) -4*p[2]]

Jf4([0.75,0.25],0.75)

function Df(p,x,M,N)
    # step 1: LD-derivatives of h1 and h1
    Dh1 = SLmin([f1(p,x) Jf1(p,x)*[M;N]],[f2(p,x) Jf2(p,x)*[M;N]])
    Dh2 = -SLmin([-f3(p,x) -Jf3(p,x)*[M;N]],[-f4(p,x) -Jf4(p,x)*[M;N]])

    # step 2: LD-derivative of f
    h1 = min(f1(p,x),f2(p,x))
    h2 = max(f3(p,x),f4(p,x))

    return SLmin([h1 Dh1],[h2 Dh2])
end


#####
# Nonsmooth newton 
function newtonSolve(p,x,M)
    rM,cM = size(M)
    
    # set up the N matrix
    N = Array{Float64,2}(undef,1,cM)
    #i = 0
    for i in 1:cM
        #i += 1
        B = [M[:,1:(i-1)] M[:,i] zeros(rM)]
        
        # Newton from here
        dk = 0.5 
        for j in 1:20

            # Forward sensitivity from here
            if i == 1
                Zk = [dk 1]
            else 
                Zk = [N[1,1:i-1] dk 1]
            end
            yk = x
            for k in 1:2
                Zk = Df(p,yk,B,Zk)
                yk = f(p,yk)
            end

            hk = dk - Zk[1,end-1]

            if abs(hk) <= 1e-8
                break
            end

            Jk = 1. - Zk[1,end]

            dk -= hk/Jk
        end

        N[1,i] = dk
    end
    return N
end




############## 
#! PART 1: make the cobweb diagrams 

#*# p2 = [0.8, 0.25]
p2 = [0.85,0.25]
x = 0:0.01:1
y = [f(p2,i) for i in x]


xc = 0.3
for i in 1:4
    yc = f1(p2,xc)
    zc = f4(p2,yc)

    dyc = Jf1(p2,xc)[3]
    dzc = Jf4(p2,yc)[3]
    dd = dzc*dyc

    xc -= (xc-zc)/(1-dd)
    println(xc-zc)
end
xC = [round(xc, digits = 16)]
for i in 1:2
    xc = round(f(p2,xc), digits = 16)
    push!(xC,xc)
end

# figure 
fig,ax = subplots(1,figsize = [5,5])
plot(x,x,color = "mediumseagreen",linewidth = 3.0,label = "\$y = x\$",ls = "--")
ax.plot(x,y,color = "darkorange",linewidth = 3.0)#, label = "y = f(p,x)")
ax.set_xlabel("\$x\$",fontsize = 18)
ax.set_ylabel("\$y\$",fontsize = 18, rotation = 0)


#CobWeb
for i in 1:length(xC)-1
    ax.plot([xC[i],xC[i],xC[i+1]],[xC[i],xC[i+1],xC[i+1]],color = "blueviolet",linewidth = 3.0,ls = ":")#,label = "2-cycle")
    if i == 1
        ax.legend(fontsize = 18)
    end   
end
ax.set_xlim([0,1])
ax.set_ylim([0,1])
fig
fig.savefig("n-cycle_2.eps")


#*# p0 = [0.75,0.25] 
p0 = [0.75,0.25]
x = 0:0.01:1
y = [f(p0,i) for i in x]


xc = 0.25
xC = [xc]
for i in 1:2
    xc = round(f(p0,xc), digits = 16)
    push!(xC,xc)
end

# figure 
fig,ax = subplots(1,figsize = [5,5])
plot(x,x,color = "mediumseagreen",linewidth = 3.0,ls = "--")#,label = "y = x")
ax.plot(x,y,color = "darkorange",linewidth = 3.0, label = "\$y = f(p,x)\$")
ax.set_xlabel("\$x\$",fontsize = 18)
ax.set_ylabel("\$y\$",fontsize = 18, rotation = 0)


#CobWeb
for i in 1:length(xC)-1
    ax.plot([xC[i],xC[i],xC[i+1]],[xC[i],xC[i+1],xC[i+1]],color = "blueviolet",linewidth = 3.0,ls = ":")#,label = "2-cycle")
    if i == 1
        ax.legend(fontsize = 18)
    end   
end
ax.set_xlim([0,1])
ax.set_ylim([0,1])

#ax.legend()
fig
fig.savefig("n-cycle_1.eps")

#*# p1 = [0.65, 0.25]
p1 = [0.65,0.25]
x = 0:0.01:1
y = [f(p1,i) for i in x]


xc = 0.3
for i in 1:4
    yc = f2(p1,xc)
    zc = f3(p1,yc)

    dyc = Jf2(p1,xc)[3]
    dzc = Jf3(p1,yc)[3]
    dd = dzc*dyc

    xc -= (xc-zc)/(1-dd)
    println(xc-zc)
end
xC = [round(xc, digits = 16)]
for i in 1:2
    xc = round(f(p1,xc), digits = 16)
    push!(xC,xc)
end

# figure 
fig,ax = subplots(1,figsize = [5,5])
plot(x,x,color = "mediumseagreen",linewidth = 3.0,ls = "--")#,label = "y = x")
ax.plot(x,y,color = "darkorange",linewidth = 3.0)#, label = "y = f(p,x)")
ax.set_xlabel("\$x\$",fontsize = 18)
ax.set_ylabel("\$y\$",fontsize = 18, rotation = 0)


#CobWeb
for i in 1:length(xC)-1
    ax.plot([xC[i],xC[i],xC[i+1]],[xC[i],xC[i+1],xC[i+1]],color = "blueviolet",linewidth = 3.0,label = "2-cycle",ls = ":")
    if i == 1
        ax.legend(fontsize = 18,loc = "lower center")
    end   
end
ax.set_xlim([0,1])
ax.set_ylim([0,1])
fig
fig.savefig("n-cycle_3.eps")



# Sensitivity Cobweb 
p0 = [0.75,0.25]
p1 = [0.65,0.25]
p2 = [0.85,0.25]
x1 = 0:0.01:1
x2 = 0:0.01:0.5
y1 = [f(p0,i) for i in x1]
y2 = [f(p1,i) for i in x2]
y3 = [f(p2,i) for i in x2]

# figure 
fig,ax = subplots(1,figsize = [6,6])
plot(x1,x1,color = "black",linewidth = 3.0)
ax.plot(x2,y2,color = "blueviolet",linewidth = 3.0,ls = "--")
ax.plot(x2,y3,color = "darkorange",linewidth = 3.0,ls = "--")
#ax.plot(x1,y1,color = "black",linewidth = 3.0)
#fig



#Cobweb calculations: 
xc = 0.25
xC1 = [xc]
for i in 1:2
    xc = round(f(p0,xc), digits = 16)
    push!(xC1,xc)
end

xc = 0.3
for i in 1:4
    yc = f2(p1,xc)
    zc = f3(p1,yc)

    dyc = Jf2(p1,xc)[3]
    dzc = Jf3(p1,yc)[3]
    dd = dzc*dyc

    xc -= (xc-zc)/(1-dd)
    println(xc-zc)
end
xC2 = [round(xc, digits = 16)]
for i in 1:2
    xc = round(f(p1,xc), digits = 16)
    push!(xC2,xc)
end

xc = 0.3
for i in 1:4
    yc = f1(p2,xc)
    zc = f4(p2,yc)

    dyc = Jf1(p2,xc)[3]
    dzc = Jf4(p2,yc)[3]
    dd = dzc*dyc

    xc -= (xc-zc)/(1-dd)
    println(xc-zc)
end
xC3 = [round(xc, digits = 16)]
for i in 1:2
    xc = round(f(p2,xc), digits = 16)
    push!(xC3,xc)
end

#CobWeb Plot
for i in 1:length(xC1)-1
    ax.plot([xC1[i],xC1[i],xC1[i+1]],[xC1[i],xC1[i+1],xC1[i+1]],color = "mediumseagreen",linewidth = 2.)
    ax.plot([xC2[i],xC2[i],xC2[i+1]],[xC2[i],xC2[i+1],xC2[i+1]],color = "blueviolet",linewidth = 2.0,ls = "--")
    ax.plot([xC3[i],xC3[i],xC3[i+1]],[xC3[i],xC3[i+1],xC3[i+1]],color = "darkorange",linewidth = 2.0,ls = "--")
end
#ax.legend()
fig

#ax.plot(x2,y2,color = "blueviolet",linewidth = 3.0,ls = "--")
#ax.plot(x2,y3,color = "darkorange",linewidth = 3.0,ls = "--")
ax.plot(x1,y1,color = "black",linewidth = 3.0)
ax.set_xlabel("x",fontsize = 16)
ax.set_ylabel("y",fontsize = 16, rotation = 0)
fig
fig.savefig("n-cycle_4.eps")








#############################################################
#############################################################
#############################################################

###! PART 2: sensitivity diagrams: 
#*# x1
p = [0.75,0.25]
Mp = [1 0;0 1]
Mm = [-1 0;0 -1]

hp = newtonSolve(p,0.25,Mp)*inv(Mp)
hm = newtonSolve(p,0.25,Mm)*inv(Mm)

pPlot = [-0.1,0,0.1] .+ 0.75
xPlot = [-hm[1]*0.1,0,hp[1]*0.1] .+ 0.25


h57 = Array{Float64,1}(undef,0)

dp = -0.1:0.001:0.1
for i in dp
    p = [0.75,0.25] + i*[1,0]#[1/sqrt(2),1/sqrt(2)]
    x = 0.15:0.000001:0.35
    y = [f(p,i) for i in x]
    z = [f(p,i) for i in y]

    h37 = x[findall(x->x<1e-5,abs.(x-z))][1]
    push!(h57,h37)
end

fig,ax = subplots(1,figsize = [8,6]*0.75)
ax.plot(dp .+ 0.75,h57,linewidth = 3,color = "black")
ax.plot(pPlot,xPlot,linewidth = 3,ls = ":",color = "red")
#ax.plot([-0.05,0,0.05],h57[1:50:end])
ax.set_ylabel("\$x_1^*\$",fontsize = 20, rotation = 0)
ax.set_xlabel("\$p_1\$",fontsize = 20, rotation = 0)
ax.set_yticks(0.22:0.025:0.32)
fig
fig.savefig("n-cycle_dir_1.eps",bbox_inches = "tight")


#*# x2
p = [0.75,0.25]
Mp = [1 0;0 1]
Mm = [-1 0;0 -1]

hp = newtonSolve(p,0.75,Mp)*inv(Mp)
hm = newtonSolve(p,0.75,Mm)*inv(Mm)

pPlot = [-0.1,0,0.1] .+ 0.75
xPlot = [-hm[1]*0.1,0,hp[1]*0.1] .+ 0.75


h57 = Array{Float64,1}(undef,0)

dp = -0.1:0.001:0.1
for i in dp
    p = [0.75,0.25] + i*[1,0]#[1/sqrt(2),1/sqrt(2)]
    x = 0.65:0.000001:0.85
    y = [f(p,i) for i in x]
    z = [f(p,i) for i in y]

    h37 = x[findall(x->x<1e-5,abs.(x-z))][1]
    push!(h57,h37)
end

fig,ax = subplots(1,figsize = [8,6]*0.75)
ax.plot(dp .+ 0.75,h57,linewidth = 3,color = "black")
ax.plot(pPlot,xPlot,linewidth = 3,ls = ":",color = "red")
#ax.plot([-0.05,0,0.05],h57[1:50:end])
ax.set_ylabel("\$x_2^*\$",fontsize = 20, rotation = 0)
ax.set_xlabel("\$p_1\$",fontsize = 20, rotation = 0)
ax.set_yticks(0.73:0.01:0.77)
fig
fig.savefig("n-cycle_dir_2.eps",bbox_inches = "tight")




