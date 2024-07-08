using PyPlot, LinearAlgebra

#### Region plot

# Region 1: 
ps111 = 0:0.1:1 
ps211 = 1 .- ps111
ps221 = [1,2]
ps121 = [0,0]


# Region 2:
ps112 = [-1,0]
ps212 = [1,1]
ps122 = [0,0]
ps222 = [0,1]


# Region 3:
ps213 = 0:0.1:1
ps113 = 1 .- ps213 
ps123 = [0,0]
ps223 = [0,1]


# Region 4:
ps114 = -1:0.001:0
ps214 = 1 ./(1 .+ps114)
ps124 = [0,0]
ps224 = [1,2]

# Region 5:
ps115 = -1:0.001:0
ps215 = 1 ./(1 .+ps114)
ps125 = [-1,0]
ps225 = [1,1]

# Region 1: 
p11 = 0:0.1:1
p211 = 1 .-p11 
p221 = 2*ones(length(p11))

# Region 2:
p12 = [-1,0]
p212 = [0,0]
p222 = [1,1]

# Region 3:
p13 = [0,1]
p213 = [0,0]
p223 = [1,0]

# Region 4:
p14 = -0.5:0.001:0
p214 = 1 ./(1 .+p14)
p224 = 2*ones(length(p14))

# Region 5:
p15 = -1.0:0.001:0
p215 = ones(length(p15))
p225 = [2*ones(500);p214]



wid = 3
### figure
fig,ax = subplots(1,figsize = [6,6])
## Region 1:
ax.fill_between(p11,p211,p221,alpha = 0.5)
ax.plot(ps111,ps211,color = "black",linewidth = wid)
ax.plot(ps121,ps221,color = "black",linewidth = wid)

## Region 2:
ax.fill_between(p12,p212,p222,alpha = 0.5)
ax.plot(ps112,ps212,color = "black",linewidth = wid)
ax.plot(ps122,ps222,color = "black",linewidth = wid)

## Region 3:
ax.fill_between(p13,p213,p223,alpha = 0.5)

## Region 4:
ax.fill_between(p14,p214,p224,alpha = 0.5)
ax.plot(ps114,ps214,color = "black",linewidth = wid)

## Region 5:
ax.fill_between(p15,p215,p225,alpha = 0.5)


ax.set_xlim([-1,1])
ax.set_ylim([0,2])

ax.set_xlabel("\$p_1\$",fontsize = 22)
ax.set_ylabel("\$p_2\$",fontsize = 22)

fig
fig.savefig("exp_1.pdf")


#########
## Orbit Plots
f(z,p) = [min(z[1],z[2]+p[1]),p[2]*z[1]]

z0 = [1,1.]
p = [0.2, 1.6]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

fig4,ax4 = subplots(2, figsize = [8,6])
ax4[1].plot(0:length(zL)-1,z1,color = "tab:blue",linewidth = wid)
ax4[1].scatter(0:length(zL)-1,z1,color = "tab:blue",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax4[2].plot(0:length(zL)-1,z2,color = "tab:blue",linewidth = wid)
ax4[2].scatter(0:length(zL)-1,z2,color = "tab:blue",marker = "o")

ax4[1].set_ylabel("\$x\$",fontsize = 22, rotation = 0)
ax4[2].set_ylabel("\$y\$",fontsize = 22, rotation = 0)
ax4[2].set_xlabel("\$k\$",fontsize = 22)

ax4[1].set_xlim([0,10])
ax4[2].set_xlim([0,10])
ax4[1].set_ylim([-0.4,1.8])
ax4[2].set_ylim([-0.4,1.8])

fig4
fig4.savefig("exp_2.eps")



z0 = [1,1.]
p = [-0.2, 0.4]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

fig4,ax4 = subplots(2, figsize = [8,6])
ax4[1].plot(0:length(zL)-1,z1,color = "tab:orange",linewidth = wid)
ax4[1].scatter(0:length(zL)-1,z1,color = "tab:orange",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax4[2].plot(0:length(zL)-1,z2,color = "tab:orange",linewidth = wid)
ax4[2].scatter(0:length(zL)-1,z2,color = "tab:orange",marker = "o")

ax4[1].set_ylabel("\$x\$",fontsize = 22, rotation = 0)
ax4[2].set_ylabel("\$y\$",fontsize = 22, rotation = 0)
ax4[2].set_xlabel("\$k\$",fontsize = 22)

ax4[1].set_xlim([0,10])
ax4[2].set_xlim([0,10])
ax4[1].set_ylim([-0.4,1.8])
ax4[2].set_ylim([-0.4,1.8])

fig4
fig4.savefig("exp_3.eps")



#### Orbit Plot 3: 
z0 = [1,1.]
p = [0.2, 1.6]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

fig5,ax5 = subplots(2,2, figsize = [10,8]*0.8)
ax5[1,1].plot(0:length(zL)-1,z1,color = "tab:blue",linewidth = wid)
ax5[1,1].scatter(0:length(zL)-1,z1,color = "tab:blue",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax5[2,1].plot(0:length(zL)-1,z2,color = "tab:blue",linewidth = wid,label = "p = (0.2,1.6)")
ax5[2,1].scatter(0:length(zL)-1,z2,color = "tab:blue",marker = "o")

ax5[1,1].set_ylabel("\$x\$",fontsize = 22, rotation = 0)
ax5[2,1].set_ylabel("\$y\$",fontsize = 22, rotation = 0)
ax5[2,1].set_xlabel("\$k\$",fontsize = 22)

ax5[1,1].set_xlim([0,10])
ax5[2,1].set_xlim([0,10])
ax5[1,1].set_ylim([-0.4,1.8])
ax5[2,1].set_ylim([-0.4,1.8])


z0 = [1,1.]
p = [-0.2, 0.4]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

ax5[1,2].plot(0:length(zL)-1,z1,color = "tab:orange",linewidth = wid)
ax5[1,2].scatter(0:length(zL)-1,z1,color = "tab:orange",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax5[2,2].plot(0:length(zL)-1,z2,color = "tab:orange",linewidth = wid,label = "p = (-0.2,0.4)")
ax5[2,2].scatter(0:length(zL)-1,z2,color = "tab:orange",marker = "o")

ax5[1,2].set_ylabel("\$x\$",fontsize = 22, rotation = 0)
ax5[2,2].set_ylabel("\$y\$",fontsize = 22, rotation = 0)
ax5[2,2].set_xlabel("\$k\$",fontsize = 22)

ax5[1,2].set_xlim([0,10])
ax5[2,2].set_xlim([0,10])
ax5[1,2].set_ylim([-0.4,1.8])
ax5[2,2].set_ylim([-0.4,1.8])

ax5[2,1].legend(fontsize = 16)
ax5[2,2].legend(fontsize = 16)
fig5
fig5.savefig("exp_4.eps")



################ Orbit plot 5
z0 = [1,1.]
p = [0.2, 1.6]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

fig6,ax6 = subplots(2, figsize = [6,6])
ax6[1].plot(0:length(zL)-1,z1,color = "tab:blue",linewidth = wid)
ax6[1].scatter(0:length(zL)-1,z1,color = "tab:blue",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax6[2].plot(0:length(zL)-1,z2,color = "tab:blue",linewidth = wid,label = "p = (0.2,1.6)")
ax6[2].scatter(0:length(zL)-1,z2,color = "tab:blue",marker = "o")


z0 = [1,1.]
p = [-0.2, 0.4]

zk = z0
zL = [zk]

for i in 1:10
    zk = f(zk,p)
    push!(zL,zk)
end

z1 = [i[1] for i in zL]
z2 = [i[2] for i in zL]

ax6[1].plot(0:length(zL)-1,z1,color = "tab:orange",linewidth = wid)
ax6[1].scatter(0:length(zL)-1,z1,color = "tab:orange",marker = "o")
#ax4[1].plot(1:length(zL),1:-0.1:0)
ax6[2].plot(0:length(zL)-1,z2,color = "tab:orange",linewidth = wid,label = "p = (-0.2,0.4)")
ax6[2].scatter(0:length(zL)-1,z2,color = "tab:orange",marker = "o")

ax6[1].set_ylabel("\$x\$",fontsize = 22, rotation = 0)
ax6[2].set_ylabel("\$y\$",fontsize = 22, rotation = 0)
ax6[2].set_xlabel("\$k\$",fontsize = 22)

ax6[1].set_xlim([0,10])
ax6[2].set_xlim([0,10])
ax6[1].set_ylim([-0.4,1.8])
ax6[2].set_ylim([-0.4,1.8])
ax6[2].legend(fontsize = 16)

fig6
fig6.savefig("exp_5.eps")