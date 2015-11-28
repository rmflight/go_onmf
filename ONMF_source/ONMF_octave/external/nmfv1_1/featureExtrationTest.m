function [testExtr,outTest,time]=featureExtrationTest(trainSet,testSet,outTrain)
% map the test/unknown data into the feature space produced by function "featureExtractionTrain".
% testSet: matrix, each column is a test/unknown sample.
% outTrain: the output of "featureExtractionTrain" function.
% testExtr: the test/unknown samples in the feature space.
% outTest: [], the other outputs, reserved for future use.
% Reference:
%  [1]\bibitem{bibm10}
%     Y. Li and A. Ngom,
%     ``Non-negative matrix and tensor factorization based classification of clinical microarray gene expression data,''
%     {\it IEEE International Conference on Bioinformatics \& Biomedicine},
%     2010, pp.438-443.
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 22, 2011

tStart=tic;
switch outTrain.feMethod
    case 'non'
        testExtr=testSet;
        outTest=[];
    case 'nmf'
        outTrain=rmfield(outTrain,'feMethod');
        testExtr=nmfnnlstest(testSet,outTrain);
        outTest=[];
    case 'sparsenmf'
        outTrain=rmfield(outTrain,'feMethod');
        testExtr=sparsenmfnnlstest(testSet,outTrain);
        outTest=[];
    case 'orthnmf'
        outTrain=rmfield(outTrain,'feMethod');
        testExtr=orthnmfruletest(testSet,outTrain);
        outTest=[];
    case 'knmf-nnls'
        AtA=outTrain.factors{1};
        trainExtr=outTrain.factors{2};
        XtS=computeKernelMatrix(trainSet,testSet,outTrain.optionTr);
        AtS=pinv(trainExtr)'*XtS;
        testExtr=kfcnnls2(AtA,AtS);
        outTest=[];
    case 'knmf-ur'
        AtA=outTrain.factors{1};
        trainExtr=outTrain.factors{2};
        XtS=computeKernelMatrix(trainSet,testSet,outTrain.optionTr);
        StS=computeKernelMatrix(testSet,testSet,outTrain.optionTr);
        AtS=pinv(trainExtr)'*XtS;
        testExtr=kurnnls(StS,AtA,AtS);
        outTest=[];    
    case 'knmf-dc'
       % the input is not testSet, it is a kernel matrix
       XtS=computeKernelMatrix(trainSet,testSet,outTrain.optionTr);
        [testExtr,~,~]=nmfnnlstest(XtS,outTrain);
        outTest=[];    
    case 'knmf-cv'
        AtA=outTrain.factors{1};
        W=outTrain.factors{2};
        trainSet=outTrain.factors{4};
        XtS=computeKernelMatrix(trainSet,testSet,outTrain.optionTr);
        AtS=W'*XtS;
        testExtr=kfcnnls2(AtA,AtS);
        outTest=[]; 
    case 'seminmf'
        outTrain=rmfield(outTrain,'feMethod');
        testExtr=seminmfnnlstest(testSet,outTrain);
        outTest=[];
    case 'convexnmf'
        outTrain=rmfield(outTrain,'feMethod');
        testExtr=convexnmfruletest(testSet,outTrain);
        outTest=[];
    case 'pca'
        testExtr=(outTrain.factors{1})'*testSet;
        outTest=[];
end
time=toc(tStart);
end