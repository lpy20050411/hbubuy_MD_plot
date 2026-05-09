#!/bin/bash

  gmx editconf -f conf.gro -o box.gro -c -d 1.2
  
  ###################
  
  gmx solvate -cp box.gro -p topol.top -o sol.gro
  
  gmx grompp -f em_steep.mdp -c sol.gro -r sol.gro -p topol.top -o ions.tpr -maxwarn 2
  
  echo SOL | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.5
  
  #########em#########
  
  gmx grompp -f em_steep.mdp -c ions.gro -r ions.gro -p topol.top -o em_steep.tpr -maxwarn 2
  
  gmx mdrun -v  -deffnm em_steep -ntmpi 1 -nt 100 -nb gpu
  
  gmx grompp -f em_cg.mdp -c em_steep.gro -r em_steep.gro -p topol.top -o em_cg.tpr -maxwarn 2
  
  gmx mdrun -v -deffnm em_cg -ntmpi 1 -nt 100 -nb gpu
  
  ########NVT#########

gmx grompp -f nvt.mdp -c em_cg.gro -r em_cg.gro -p topol.top -o nvt.tpr -maxwarn 2

gmx mdrun -v  -deffnm nvt -ntmpi 1 -nt 100 -nb gpu

gmx mdrun -v  -deffnm nvt -ntmpi 1 -nt 24 -nb gpu -pme gpu -bonded gpu


#########NPT########

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 2 

gmx mdrun -v  -deffnm npt -ntmpi 1 -nt 100 -nb gpu

gmx mdrun -v  -deffnm npt -ntmpi 1 -nt 24 -nb gpu -pme gpu -bonded gpu


########MD##########

gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 2

gmx mdrun -v  -deffnm md -ntmpi 1 -nt 100 -nb gpu

gmx mdrun -v  -deffnm md -ntmpi 1 -nt 24 -nb gpu -pme gpu -bonded gpu



####analysis#########

echo "q" | gmx make_ndx -f md.tpr -o sys.ndx

echo -e "MOL \n 0" | gmx trjconv -f md.xtc -o md_1.xtc -s md.tpr -n sys.ndx -pbc mol -center -ur compact -skip 10

echo -e "MOL \n 0" | gmx trjconv -f md_1.xtc -o md_2.xtc -s md.tpr -pbc nojump

echo -e "21 \n 21" | gmx rms -s md.tpr -n sys.ndx -f md.xtc -o RMSD_Ca.xvg   

echo -e "21 " | gmx gyrate -s md.tpr -f md.xtc -o Gyrate.xvg                   backbone

echo -e "21 " | gmx rmsf -s md.tpr -n sys.ndx -f md.xtc -o RMSF.xvg            C-alpha



steps=(0.02 0.01 0.005 0.002)

for i in ${!steps[@]}
do
  step=${steps[$i]}
  echo "=== Run $i with emstep=$step ==="

  # 生成新的 mdp
  sed "s/^emstep.*/emstep = $step/" em_steep.mdp > em_steep_$i.mdp

或者

for i in {1..5}
do


断跑继续（不新生成xtc，不延长时间）
nohup mpirun -np 1 gmx_mpi mdrun \
  -deffnm md \
  -cpi md.cpt \
  -append \
  -v

延长模拟时间convert-tpr

 第一步：延长 tpr
  gmx convert-tpr \
  -s md.tpr \
  -extend 50000 \
  -o md_extend.tpr
 第二步：继续跑（重点）
  gmx_mpi mdrun \
  -s md_extend.tpr \
  -deffnm md \
  -cpi md.cpt \
  -append \
  -v
  -np 1

实测4090不用-update gpu