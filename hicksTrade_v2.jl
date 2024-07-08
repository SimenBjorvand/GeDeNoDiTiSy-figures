using PyPlot, LinearAlgebra

# Model: 
function f(x,p)
    y,yp = x
    a,b,c,d = p

    v = min(b,max(a*(y-yp),-d))
    yn = c*y + v
    ypn = y

    return [yn,ypn],v
end

# Reverse model
function frev(x,p)
    y,yp = x
    a,b,c,d = p

    yn = yp
    ypn = (a+c)/a*yp - y/a
    return [yn,ypn]
end

####################################
####################################

# LD-derivatives
function lmin(a::Array{Float64,1},b::Array{Float64,1})
    n = length(a)
    for i in 1:n
        if a[i] < b[i]
            return a
        elseif b[i] < a[i]
            return b
        end
    end
    return a  
end

function slmin(a::Array{Float64,1},b::Array{Float64,1})
    m = lmin(a,b)
    return m[2:end]
end
function slmax(a::Array{Float64,1},b::Array{Float64,1})
    return -slmin(-a,-b)
end


function fLD(x,p,X,M)
    y,yp = x
    a,b,c,d = p
    ~,q = size(X)

    Y = X[1,:]
    Yp = X[2,:]
    Ma = M[1,:]
    Mb = M[2,:]
    Mc = M[3,:]
    Md = M[4,:]

    Xn = zeros(2,q)

    U1 = [a*(y-yp);(y-yp)*Ma+a*(Y-Yp)]
    U2 = [-d;-Md]
    U = slmax(U1,U2)

    V1 = [b;Mb]
    V2 = [max(a*(y-yp),-d);U]
    V = slmin(V1,V2)

    Xn[1,:] = c*Y + y*Mc + V
    Xn[2,:] = Y

    return Xn,V
end

# Consider only the case where c < 1 and set b = d as a new parameter
function fLDr(x,pr,X,M)
    y,yp = x
    a,b,c = pr
    ~,q = size(X)

    Y = X[1,:]
    Yp = X[2,:]
    Ma = M[1,:]
    Mb = M[2,:]
    Mc = M[3,:]

    Xn = zeros(2,q)

    U1 = [a*(y-yp);(y-yp)*Ma+a*(Y-Yp)]
    U2 = [-b;-Mb]
    U = slmax(U1,U2)

    V1 = [b;Mb]
    V2 = [max(a*(y-yp),-b);U]
    V = slmin(V1,V2)

    Xn[1,:] = c*Y + y*Mc + V
    Xn[2,:] = Y

    return Xn,V
end


#############
# Figure:
p = [0.9,0.5,0.5,0.5]
pr = [0.9,0.5,0.5]

x0 = [0.75,0.75]
M = I(3)
Xk = zeros(2,3)

xL = [x0]
zL = [0.0]
X1L = [zeros(3)]
ZL = [zeros(3)]

xk = x0

xLimit = 20
for i in 1:xLimit
    Xk,Zk = fLDr(xk,pr,Xk,M)
    xk,zk = f(xk,p)

    push!(xL,xk)
    push!(zL,zk)
    push!(X1L,Xk[1,:])
    push!(ZL,Zk)
end

x1 = [i[1] for i in xL]
zL
x1a = [i[1] for i in X1L]
x1b = [i[2] for i in X1L]
x1c = [i[3] for i in X1L]
za = [i[1] for i in ZL]
zb = [i[2] for i in ZL]
zc = [i[3] for i in ZL]

fig,ax = subplots(3,2,figsize = [12,8])
ax[1,1].plot(0:length(x1)-1,x1)
ax[1,2].plot(0:length(zL)-1,zL)
ax[1,1].scatter(0:length(x1)-1,x1,marker = "o")
ax[1,2].scatter(0:length(zL)-1,zL,marker = "o")

ax[2,1].plot(0:length(x1)-1,x1a)
ax[3,1].plot(0:length(x1)-1,x1b,label = "a < 1")
ax[2,2].plot(0:length(zL)-1,za)
ax[3,2].plot(0:length(zL)-1,zb)
ax[2,1].scatter(0:length(x1)-1,x1a,marker = "o")
ax[3,1].scatter(0:length(x1)-1,x1b,marker = "o")
ax[2,2].scatter(0:length(zL)-1,za,marker = "o")
ax[3,2].scatter(0:length(zL)-1,zb,marker = "o")



# 2
p = [1.1,0.5,0.5,0.5]
pr = [1.1,0.5,0.5]

x0 = [0.75,0.75]
M = I(3)
Xk = zeros(2,3)

xL = [x0]
zL = [0.0]
X1L = [zeros(3)]
ZL = [zeros(3)]

xk = x0

xLimit = 20
for i in 1:xLimit
    Xk,Zk = fLDr(xk,pr,Xk,M)
    xk,zk = f(xk,p)

    push!(xL,xk)
    push!(zL,zk)
    push!(X1L,Xk[1,:])
    push!(ZL,Zk)
end

x1 = [i[1] for i in xL]
zL
x1a = [i[1] for i in X1L]
x1b = [i[2] for i in X1L]
x1c = [i[3] for i in X1L]
za = [i[1] for i in ZL]
zb = [i[2] for i in ZL]
zc = [i[3] for i in ZL]


ax[1,1].plot(0:length(x1)-1,x1)
ax[1,2].plot(0:length(zL)-1,zL)
ax[1,1].scatter(0:length(x1)-1,x1,marker = "o")
ax[1,2].scatter(0:length(zL)-1,zL,marker = "o")

ax[2,1].plot(0:length(x1)-1,x1a)
ax[3,1].plot(0:length(x1)-1,x1b,label = "a > 1")
ax[2,2].plot(0:length(zL)-1,za)
ax[3,2].plot(0:length(zL)-1,zb)
ax[2,1].scatter(0:length(x1)-1,x1a,marker = "o")
ax[3,1].scatter(0:length(x1)-1,x1b,marker = "o")
ax[2,2].scatter(0:length(zL)-1,za,marker = "o")
ax[3,2].scatter(0:length(zL)-1,zb,marker = "o")

ax[1,1].set_ylim([-1.0,1])
ax[1,1].set_xlim([0,xLimit])
ax[1,2].set_ylim([-1,1])
ax[1,2].set_xlim([0,xLimit])


ax[2,1].set_xlim([0,xLimit])
ax[2,2].set_xlim([0,xLimit])
ax[3,1].set_xlim([0,xLimit])
ax[3,2].set_xlim([0,xLimit])
ax[2,1].set_ylim([-3.4,3.4])
ax[2,2].set_ylim([-3.4,3.4])
ax[3,1].set_ylim([-3.4,3.4])
ax[3,2].set_ylim([-3.4,3.4])


ax[1,1].set_ylabel("\$y\$",fontsize = 16, rotation = 0)
ax[1,2].set_ylabel("\$v\$",fontsize = 16, rotation = 0)
ax[2,1].set_ylabel("\$S_{y_a}\$",fontsize = 16, rotation = 0)
ax[2,2].set_ylabel("\$S_{v_a}\$",fontsize = 16, rotation = 0)
ax[3,1].set_ylabel("\$S_{y_b}\$",fontsize = 16, rotation = 0)
ax[3,2].set_ylabel("\$S_{v_b}\$",fontsize = 16, rotation = 0)


ax[3,1].set_xlabel("\$k\$",fontsize = 16)
ax[3,2].set_xlabel("\$k\$",fontsize = 16)

ax[1,1].set_xticks(0:4:20)
ax[2,1].set_xticks(0:4:20)
ax[3,1].set_xticks(0:4:20)

ax[1,2].set_xticks(0:4:20)
ax[2,2].set_xticks(0:4:20)
ax[3,2].set_xticks(0:4:20)

ax[3,1].legend()


fig
fig.savefig("hicks_sim.eps")