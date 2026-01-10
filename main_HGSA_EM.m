%% Information
% Chauhan, Dikshit. "Restart mechanism-based multilevel gravitational search algorithm for global optimization and image segmentation." 
% Engineering Applications of Artificial Intelligence 163 (2026): 112904.
% DOI: https://doi.org/10.1016/j.engappai.2025.112904

clear all;clc
warning off

fhd = str2func('cec17_func');
% D=30;N=2*D;max_fe=4000*D;
% layers=[4 10 18 28];
D=50;N=2*D;max_fe=4000*D;
layers=[2 4 8 18 32 40];
JJ=[1,3:30];
for i = JJ
    func_num=i
    [Fbest_HGSA_EM,BestChart_HGSA_EM]=MGSA_EM(fhd,layers,D,max_fe,func_num);Fbest_HGSA_EM;
   %% figure
    semilogy(BestChart_HGSA_EM,'Color','r','LineWidth',2.5)
    xlabel('Iteration numbers','Fontname','Times New Roman','fontsize',12','FontWeight','bold');
    ylabel('Best values so far','Fontname','Times New Roman','fontsize',12,'FontWeight','bold');
    axis tight
    grid on
    set(gca,'fontweight','bold','fontsize',12,'fontname','Times New Roman')
    box on
    legend('\fontsize{10}\bf MGSA-EM')
end
