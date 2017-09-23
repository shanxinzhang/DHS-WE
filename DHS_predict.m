function issuccess=DHS_predict(inputfile,outputfile)
% Predictiion DHSs sequences from file.
%
%   DHS_predict(inputfile,outputfile) predict sequences from inputfile, and write results to outputfile

%  issucess=DHS_predict(inputfile,outputfile), predict sequences from
%  inputfile, and write results to outputfile, if success, return 1, else
%  return 0

%   DHS_predict(inputfile) predict sequences from inputfile, and write
%   results to 'pred_result.txt' file.

%   Notes:
%
%   The inputfile is a fasta format file.
%
%   Examples:
%
% issucess=DHS_predict('test.fasta','test_out.txt')
%   

% Shanxin Zhang
% shanxinzhang@jiangnan.edu.cn

if nargin<1
    issuccess=0;
    disp('No input files found');
    return;
elseif nargin<2
    outputfile='pred_result.txt';
end
%please replace the direction of python in your system
command=['C:\Python27\python.exe GenerateFeatures.py ', inputfile];
system(command);

[head,~]=fastaread(inputfile);
labels=ones(length(head),1);
    for i=2:6
        kmer_data{i-1}=importdata(['DHSs_kmer_',num2str(i),'.txt']);
        reckmer_data{i-1}=importdata(['DHSs_reckmer_',num2str(i),'.txt']);
    end
        for i=3:6
         m_data=importdata(['DHSs_MismatchProfile_',num2str(i),'_1.txt']);
         mismatch_data{i-2}=m_data./repmat(sum(m_data,2),1,size(m_data,2));
        end
        l=[2,1,6,6];
        w=[0.2,0.1,0.1,0.1];
        for i=1:4
             if i==1
                pseknc_data{i}=importdata(['DHSs_PC_PseDNC_',num2str(l(i)),'_',num2str(w(i)),'.txt']);
            elseif i==2
                pseknc_data{i}=importdata(['DHSs_SC_PseDNC_',num2str(l(i)),'_',num2str(w(i)),'.txt']);       
            elseif i==3
                pseknc_data{i}=importdata(['DHSs_PC_PseTNC_',num2str(l(i)),'_',num2str(w(i)),'.txt']);
            else
                pseknc_data{i}=importdata(['DHSs_SC_PseTNC_',num2str(l(i)),'_',num2str(w(i)),'.txt']);
             end
        end
        w=[0.5,0.7];
        for i=2:3
            j=w(i-1);
            s_data=importdata(['DHSs_SubsequenceProfile_',num2str(i),'_',num2str(j),'.txt']);
            subseq_data{i-1}=s_data./repmat(sum(s_data,2),1,size(s_data,2));
        end
            kacc_data=importdata('DHSs_dacc_1.txt');

           cur_mat_file=dir([pwd,'\models\','*.mat']);
           for i=1:length(cur_mat_file)
               load([pwd,'\models\',cur_mat_file(i).name]);
           end
             final_models_14={final_models_14_1;final_models_14_2;final_models_14_3;final_models_14_4;final_models_14_5;final_models_14_6;...
                 final_models_14_7;final_models_14_8;final_models_14_9;final_models_14_10;final_models_14_11;final_models_14_12;...
                 final_models_14_13;final_models_14_14;final_models_14_15};

           for i=1:15
              [~,~,pred_result]=libsvmpredict(labels,kmer_data{1},final_models_1{i,1},'-b 1 -q');pred_prob(:,1)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,kmer_data{2},final_models_2{i,1},'-b 1 -q');pred_prob(:,2)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,kmer_data{3},final_models_3{i,1},'-b 1 -q');pred_prob(:,3)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,kmer_data{4},final_models_4{i,1},'-b 1 -q');pred_prob(:,4)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,kmer_data{5},final_models_5{i,1},'-b 1 -q');pred_prob(:,5)=pred_result(:,1);
              
              [~,~,pred_result]=libsvmpredict(labels,reckmer_data{1},final_models_6{i,1},'-b 1 -q');pred_prob(:,6)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,reckmer_data{2},final_models_7{i,1},'-b 1 -q');pred_prob(:,7)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,reckmer_data{3},final_models_8{i,1},'-b 1 -q');pred_prob(:,8)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,reckmer_data{4},final_models_9{i,1},'-b 1 -q');pred_prob(:,9)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,reckmer_data{5},final_models_10{i,1},'-b 1 -q');pred_prob(:,10)=pred_result(:,1);      
              
              [~,~,pred_result]=libsvmpredict(labels,mismatch_data{1},final_models_11{i,1},'-b 1 -q');pred_prob(:,11)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,mismatch_data{2},final_models_12{i,1},'-b 1 -q');pred_prob(:,12)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,mismatch_data{3},final_models_13{i,1},'-b 1 -q');pred_prob(:,13)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,mismatch_data{4},final_models_14{i,1},'-b 1 -q');pred_prob(:,14)=pred_result(:,1);  

              [~,~,pred_result]=libsvmpredict(labels,subseq_data{1},final_models_15{i,1},'-b 1 -q');pred_prob(:,15)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,subseq_data{2},final_models_16{i,1},'-b 1 -q');pred_prob(:,16)=pred_result(:,1);
              
              [~,~,pred_result]=libsvmpredict(labels,kacc_data,final_models_17{i,1},'-b 1 -q');pred_prob(:,17)=pred_result(:,1);     
              
              [~,~,pred_result]=libsvmpredict(labels,pseknc_data{1},final_models_18{i,1},'-b 1 -q');pred_prob(:,18)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,pseknc_data{2},final_models_19{i,1},'-b 1 -q');pred_prob(:,19)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,pseknc_data{3},final_models_20{i,1},'-b 1 -q');pred_prob(:,20)=pred_result(:,1);
              [~,~,pred_result]=libsvmpredict(labels,pseknc_data{4},final_models_21{i,1},'-b 1 -q');pred_prob(:,21)=pred_result(:,1);
                            
              final_pred_result(:,i)=pred_prob*GA_final_w{i}';
           end
           pred_label=sign(mean(sign(final_pred_result-0.5),2));
           delete('DHS*');
           %outputfile='test_out.txt';
           fid=fopen(outputfile,'w');
        for i=1:length(head)
            if pred_label(i)==1
                fprintf(fid,'%s is predicted as DHSs.\n',head{i});
            else
                fprintf(fid,'%s is predicted as Non DHSs.\n',head{i});
            end
        end
        fclose(fid);
        clear all;
        issuccess=1;
end