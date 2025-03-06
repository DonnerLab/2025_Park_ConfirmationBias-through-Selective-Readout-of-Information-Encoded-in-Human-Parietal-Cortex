function [groi,gname, gid] = glasser_group(indx)
%% 22 groups of all Glasser ROIs

%%%%% DO NOT CHANGE THE ORDER OF THE PARCELS %%%%%%%


areas_grp = {};
areas_name = {};
areas_id = {};

areas_name{1} = 'V1';
areas_id{1} = 'V1';
areas_grp{1} = {'V1'};

areas_name{2} = 'Early Visual';
areas_id{2} = 'EV';
areas_grp{2} = {'V2','V3','V4'};

areas_name{3} = 'Dorsal Stream V';
areas_id{3} = 'DSV';
areas_grp{3} = {'V6','V3A','V7','IPS1','V3B','V6A'};

areas_name{4} = 'Ventral Stream V';
areas_id{4} = 'VSV';
areas_grp{4} = {'V8','FFC','PIT','VMV1','VMV3','VMV2','VVC'};

areas_name{5} = 'MT+Complex V';
areas_id{5} = 'MTCV';
areas_grp{5} = {'MST','LO1','LO2','MT','PH','V4t','FST','V3CD','LO3'};

areas_name{6} = 'Somatosensory+Motor';
areas_id{6} = 'SSM';
areas_grp{6} = {'4','3b','1','2','3a'};

areas_name{7} = 'Paracentral+Mid Cingulate';
areas_id{7} = 'PC-mCin';
areas_grp{7} = {'5m','5mv','5L','24dd','24dv','SCEF','6ma','6mp'};

areas_name{8} = 'Premotor';
areas_id{8} = 'PrM';
areas_grp{8} = {'FEF','PEF','6d','6v','6r','6a','55b'};

areas_name{9} = 'Posterior Opercular';
areas_id{9} = 'poOper';
areas_grp{9} = {'43','OP4','OP1','OP2-3','FOP1','PFcm'};

areas_name{10} = 'Early Auditory';
areas_id{10} = 'EA';
areas_grp{10} = {'A1','52','RI','PBelt','MBelt','LBelt'};

areas_name{11} = 'Auditory Association';
areas_id{11} = 'AAss';
areas_grp{11} = {'TA2','STGa','A5','STSda','STSdp','STSvp','A4','STSva'};

areas_name{12} = 'Insular+Frontal Opercular'; 
areas_id{12} = 'InFOper';
areas_grp{12} = {'PoI2','FOP4','MI','Pir','AVI','AAIC','FOP3','FOP2','PoI1','Ig','FOP5','PI'};

areas_name{13} = 'Medial Temporal';
areas_id{13} = 'MdTemp';
areas_grp{13} = {'EC','PreS','H','PeEc','PHA1','PHA3','TF','PHA2'};

areas_name{14} = 'Lateral Temporal';
areas_id{14} = 'LtTemp';
areas_grp{14} = {'TGd','TE1a','TE1p','TE2a','TE2p','PHT','TGv','TE1m'};

areas_name{15} = 'Temporo-Parietal-Occipital J';
areas_id{15} = 'TempParOcc';
areas_grp{15} = {'PSL','STV','TPOJ1','TPOJ2','TPOJ3'};

areas_name{16} = 'Sup. Parietal';
areas_id{16} = 'SupPar';
areas_grp{16} = {'7Pm','7AL','7Am','7PL','7PC','LIPv','VIP','MIP','LIPd','AIP'};

areas_name{17} = 'Inferior Parietal';
areas_id{17} = 'InfPar';
areas_grp{17} = {'PFt','PGp','IP2','IP1','IP0','PFop','PF','PFm','PGi','PGs'};

areas_name{18} = 'Post. Cingulate';
areas_id{18} = 'pCin';
areas_grp{18} = {'RSC','POS2','PCV','7m','POS1','23c','23d','v23ab','d23ab','31pv','ProS','DVT','31pd','31a'};

areas_name{19} = 'Anterior Cingulate+Medial Prefrontal';
areas_id{19} = 'aCinMdPFC';
areas_grp{19} = {'p24pr','33pr','a24pr','p32pr','a24','d32','8BM','p32','10r','9m','10v','25','s32','a32pr','p24'};

areas_name{20} = 'Orbital+Polar Frontal'; 
areas_id{20} = 'OrbPlrFC';
areas_grp{20} = {'47m','10d','a10p','10pp','11l','13l','OFC','47s','p10p','pOFC'};

areas_name{21} = 'Inferior Frontal';
areas_id{21} = 'InfFC';
areas_grp{21} = {'44','45','47l','a47r','IFJa','IFJp','IFSp','IFSa','p47r'};

areas_name{22} = 'DL Prefrontal';
areas_id{22} = 'DLPFC';
areas_grp{22} = {'SFL','8Av','8Ad','8BL','9p','8C','p9-46v','46','a9-46v','9-46d','9a','i6-8','s6-8'};

areas_name{23} = 'Eye Fields';
areas_id{23} = 'EYE';
areas_grp{23} = {'FEF','PEF'};

areas_name{24} = 'PMv';
areas_id{24} = 'PMv';
areas_grp{24} = {'6v','6r'};

areas_name{25} = 'PMd';
areas_id{25} = 'PMd';
areas_grp{25} = {'6d','6a'};

areas_name{26} = 'M1 (hand)';
areas_id{26} = 'M1';
areas_grp{26} = {'JWG_M1'};

groi = areas_grp{indx};
gname = areas_name{indx};
gid = areas_id{indx};
return
% sum(cellfun(@(x) length(x), areas_grp))

% groups of interest
% indx = [1 2 3 5 4 16 8 13 18 22 26];
% sum(cellfun(@(x) length(x), areas_grp(indx)))
