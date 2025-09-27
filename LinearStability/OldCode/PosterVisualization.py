fig_poster, ax_poster = plt.subplots(1, 1, figsize = (4, 7),
                                subplot_kw = dict(projection = "polar"))
ax_poster.grid(False)
ax_poster.pcolormesh(phiVisCheb, rVisCheb, psiCheb.real, cmap = "RdBu_r", vmin = -1, vmax = 1)
#ax_poster.set_title(f"$k=$ {kphi}; $m=$ {kz}")
ax_poster.grid(True)

fig_poster.suptitle(f"Re [$\hat{{\psi}}(r)$ exp$(ik\phi)$];\n fastest-growing mode with\n $k=$ {kphi} and $m=$ {kz}\n\n")
fig_poster.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), cmap = "RdBu_r"),
                     ax = ax_poster, orientation = "horizontal") #, shrink = 0.75)
plt.show()
fig_poster.savefig(f"poster_streamfunc_k{kphi}_m{kz}.pdf")
plt.close(fig_poster)

xx, yy = rVisCheb * np.cos(phiVisCheb), rVisCheb * np.sin(phiVisCheb)

fig, axs = plt.subplots(1, 2, figsize = (9, 5),
                                subplot_kw = dict(projection = "3d"))

axs[0].plot_surface(xx, yy, psiCheb.real,
                            rstride = 2, cstride = 2, cmap = "RdYlBu_r",
                            vmin = -1, vmax = 1, alpha = 0.8)
axs[1].plot_surface(xx, yy, psiCheb.imag,
                            rstride = 2, cstride = 2, cmap = "RdYlBu_r",
                            vmin = -1, vmax = 1, alpha = 0.8)

for i in range(2):
    axs[i].set_xlabel("x")
    axs[i].set_ylabel("y")
    axs[i].set_xlim(-paramsCheb.Lr, paramsCheb.Lr)
    axs[i].set_ylim(-paramsCheb.Lr, paramsCheb.Lr)

axs[0].set_zlabel(f"Re[$\hat{{\psi}}(r)$]; Cheb solver")
axs[1].set_zlabel(f"Im[$\hat{{\psi}}(r)$]; Cheb solver")

fig.subplots_adjust(hspace = 0.5, wspace = 0.1)
fig.suptitle(f"Components of mode-{jj} eigen-streamfunction in $r\phi$-plane for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                    cmap = "RdYlBu_r"),
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.75)
plt.show()
fig.savefig(f"streamfunc_surface_k{kphi}_m{kz}_mode{jj}.png")
plt.close(fig)

