import numpy as np 
import matplotlib.pyplot as plt

denssize = np.loadtxt('sizedense.dat')
densdens = np.loadtxt('Desitydense.dat')
sparsesize = np.loadtxt('sizesparse.dat')
sparsedens = np.loadtxt('Desitysparse.dat')
lapackdens = np.loadtxt('DesityLAPACK.dat')
lapacksize = np.loadtxt('sizeLAPACK.dat')

fds = np.polyfit(x = np.log10(denssize[-50:-1,0]),y = np.log10(denssize[-50:-1,1]), deg=1)
fss = np.polyfit(x = np.log10(sparsesize[-20:-1,0]),y = np.log10(sparsesize[-20:-1,1]), deg=1)
fls = np.polyfit(x = np.log10(lapacksize[-50:-1,0]),y = np.log10(lapacksize[-50:-1,1]), deg=1)
FDS = np.poly1d(fds)
FSS = np.poly1d(fss)
FLS = np.poly1d(fls)
print(F'dense sizem {fds}')
print(F"sparse sizem {fss}")
print(F'lapack sizem {fls}')

fdt = np.polyfit(x = np.log10(denssize[-50:-1,0]),y = np.log10(denssize[-50:-1,2]), deg=1)
fst = np.polyfit(x = np.log10(sparsesize[-20:-1,0]),y = np.log10(sparsesize[-20:-1,2]), deg=1)
flt= np.polyfit(x = np.log10(lapacksize[-50:-1,0]),y = np.log10(lapacksize[-50:-1,2]), deg=1)
FDT = np.poly1d(fdt)
FST = np.poly1d(fst)
FLT = np.poly1d(flt)
print(F"dense sizet{fdt}")
print(F'sparse sizet{fst}')
print(f'lapack sizet{flt}')


'''
plt.plot(denssize[:-1,0],np.log10(denssize[:-1,1]), label = 'dense-Nmult')
plt.plot(sparsesize[:-1,0],np.log10(sparsesize[:-1,1]), label = 'sparce-Nmult')
plt.plot(lapacksize[:-1,0],np.log10(lapacksize[:-1,1]), label = 'LAPACK-Nmult')
plt.plot(denssize[:-1,0],FDS(denssize[:-1,0]))
plt.plot(sparsesize[:-1,0],FSS(sparsesize[:-1,0]))
plt.plot(lapacksize[:-1,0],FLS(lapacksize[:-1,0]))
plt.xlabel('matrix size')
plt.ylabel('ln(FLOPs)')
plt.legend()
plt.show()


plt.plot(denssize[:-1,0],np.log10(denssize[:-1,2]), label = 'dense-time')
plt.plot(sparsesize[:-1,0],np.log10(sparsesize[:-1,2]), label = 'sparce-time')
plt.plot(lapacksize[:-1,0],np.log10(lapacksize[:-1,2]), label = 'LAPACK-time')
plt.plot(denssize[:-1,0],FDT(denssize[:-1,0]))
plt.plot(sparsesize[:-1,0],FST(sparsesize[:-1,0]))
plt.plot(lapacksize[:-1,0],FLT(lapacksize[:-1,0]))
plt.xlabel('matrix size')
plt.ylabel('ln(time)')
plt.legend()
plt.show()
'''


# Create a figure with 2 subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Adjust `figsize` as needed

# FLOPs plot (left subplot)
axes[0].plot(np.log10(denssize[:-1,0]), np.log10(denssize[:-1,1]), label='dense-FLOPs')
axes[0].plot(np.log10(sparsesize[:-1,0]), np.log10(sparsesize[:-1,1]), label='sparse-FLOPs')
axes[0].plot(np.log10(lapacksize[:-1,0]), np.log10(lapacksize[:-1,1]), label='LAPACK-FLOPs')

axes[0].plot(np.log10(denssize[:-1,0]), FDS(np.log10(denssize[:-1,0])), color='C0', linestyle='--')  # Match dense-Nmult color
axes[0].plot(np.log10(sparsesize[:-1,0]), FSS(np.log10(sparsesize[:-1,0])), color='C1', linestyle='--')  # Match sparse-Nmult color
axes[0].plot(np.log10(lapacksize[:-1,0]), FLS(np.log10(lapacksize[:-1,0])), color='C2', linestyle='--')  # Match LAPACK-Nmult color

axes[0].set_xlabel('log(Matrix size)')
axes[0].set_ylabel('log(FLOPs)')
axes[0].legend()
axes[0].set_title('FLOPs vs Matrix Size')

# Time plot (right subplot)
axes[1].plot(np.log10(denssize[:-1,0]), np.log10(denssize[:-1,2]), label='dense-time')
axes[1].plot(np.log10(sparsesize[:-1,0]), np.log10(sparsesize[:-1,2]), label='sparse-time')
axes[1].plot(np.log10(lapacksize[:-1,0]), np.log10(lapacksize[:-1,2]), label='LAPACK-time')

axes[1].plot(np.log10(denssize[:-1,0]), FDT(np.log10(denssize[:-1,0])), color='C0', linestyle='--')  # Match dense-time color
axes[1].plot(np.log10(sparsesize[:-1,0]), FST(np.log10(sparsesize[:-1,0])), color='C1', linestyle='--')  # Match sparse-time color
axes[1].plot(np.log10(lapacksize[:-1,0]), FLT(np.log10(lapacksize[:-1,0])), color='C2', linestyle='--')  # Match LAPACK-time color

axes[1].set_xlabel('log(Matrix size)')
axes[1].set_ylabel('log(Time)')
axes[1].legend()
axes[1].set_title('Time vs Matrix Size')

# Adjust layout for better spacing
fig.tight_layout()

# Show the combined plot
plt.show()



Gds = np.polyfit(x = np.log10(densdens[-20:-1,0]),y = np.log10(densdens[-20:-1,1]), deg=1)
Gss = np.polyfit(x = np.log10(sparsedens[-20:-1,0]),y = np.log10(sparsedens[-20:-1,1]), deg=1)
Gls = np.polyfit(x = np.log10(lapackdens[-20:-1,0]),y = np.log10(lapackdens[-20:-1,1]), deg=1)
GDS = np.poly1d(Gds)
GSS = np.poly1d(Gss)
GLS = np.poly1d(Gls)
print(f"dense densen {Gds}")
print(f"sparse densen {Gss}")
print(f'l densen {Gls}')

Gdt = np.polyfit(x = np.log10(densdens[-20:-1,0]),y = np.log10(densdens[-20:-1,2]), deg=1)
Gst = np.polyfit(x = np.log10(sparsedens[-20:-1,0]),y = np.log10(sparsedens[-20:-1,2]), deg=1)
Glt= np.polyfit(x = np.log10(lapackdens[-20:-1,0]),y = np.log10(lapackdens[-20:-1,2]), deg=1)
GDT = np.poly1d(Gdt)
GST = np.poly1d(Gst)
GLT = np.poly1d(Glt)
print(f'dense denst {Gdt}')
print(f'sparse denset {Gst}')
print(f"l denst {Glt}")



import matplotlib.pyplot as plt
import numpy as np

# Create a figure with 2 subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Adjust `figsize` as needed

# FLOPs plot (left subplot)
axes[0].plot(np.log10(densdens[:-1,0]), np.log10(densdens[:-1,1]), label='dense-FLOPs')
axes[0].plot(np.log10(sparsedens[:-1,0]), np.log10(sparsedens[:-1,1]), label='sparse-FLOPs')
axes[0].plot(np.log10(lapackdens[:-1,0]), np.log10(lapackdens[:-1,1]), label='LAPACK-FLOPs')

axes[0].plot(np.log10(densdens[:-1,0]), GDS(np.log10(densdens[:-1,0])), color='C0', linestyle='--')  # Match dense-Nmult color
axes[0].plot(np.log10(sparsedens[:-1,0]), GSS(np.log10(sparsedens[:-1,0])), color='C1', linestyle='--')  # Match sparse-Nmult color
axes[0].plot(np.log10(lapackdens[:-1,0]), GLS(np.log10(lapackdens[:-1,0])), color='C2', linestyle='--')  # Match LAPACK-Nmult color

axes[0].set_xlabel('log(%Filled)')
axes[0].set_ylabel('log(FLOPs)')
axes[0].legend()
axes[0].set_title('FLOPs vs %Filled')

# Time plot (right subplot)
axes[1].plot(np.log10(densdens[:-1,0]), np.log10(densdens[:-1,2]), label='dense-time')
axes[1].plot(np.log10(sparsedens[:-1,0]), np.log10(sparsedens[:-1,2]), label='sparse-time')
axes[1].plot(np.log10(lapackdens[:-1,0]), np.log10(lapackdens[:-1,2]), label='LAPACK-time')

axes[1].plot(np.log10(densdens[:-1,0]), GDT(np.log10(densdens[:-1,0])), color='C0', linestyle='--')  # Match dense-time color
axes[1].plot(np.log10(sparsedens[:-1,0]), GST(np.log10(sparsedens[:-1,0])), color='C1', linestyle='--')  # Match sparse-time color
axes[1].plot(np.log10(lapackdens[:-1,0]), GLT(np.log10(lapackdens[:-1,0])), color='C2', linestyle='--')  # Match LAPACK-time color

axes[1].set_xlabel('log(%Filled)')
axes[1].set_ylabel('log(Time)')
axes[1].legend()
axes[1].set_title('Time vs %Filled')

# Adjust layout for better spacing
fig.tight_layout()

# Show the combined plot
plt.show()


 