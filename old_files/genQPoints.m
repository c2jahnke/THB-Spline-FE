function QPoints = genQPoints(obj,Points,lvl)
% within ImpointCurvePlot generate Points in higher basis;
        sQP = obj.levelBas{lvl}.n-obj.levelBas{lvl}.p+1;
        QPoints = zeros(sQP,2);
        QPoints(:,1) = [obj.levelBas{lvl}.a:(obj.levelBas{lvl}.b-obj.levelBas{lvl}.a)/(sQP-1):obj.levelBas{lvl}.b]'
        cnt = 1;
        for k = 1 : sQP
            if(QPoints(k,1) == Points(cnt,1))
                QPoints(k,2) = Points(cnt,2);
                cnt = cnt +1;
            end
        end
end