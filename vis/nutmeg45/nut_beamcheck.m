% function nut_beamcheck

global nuts beam;
b = squeeze(nuts.Lp(:,1,:))*beam.s_th + squeeze(nuts.Lp(:,2,:))*beam.s_ph;
figure; plot(nuts.meg.latency,b); axis tight;