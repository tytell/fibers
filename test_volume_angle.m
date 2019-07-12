function test_volume_angle

sz = 33;
a = [0.1; 0.3; 0.8];
b = [0.7; 0.1; 0.1];
c = [0; 1; 0];

A = [a b c];
A = gramschmidt(A);

a = A(:,1);
b = A(:,2);
c = A(:,3);

ctr1 = [-0.3 0.03 0.4]*sz;
ctr2 = [0.33 -0.1 0.4]*sz;
d1 = [0.15 0.15 0.15]*sz;
d2 = [0.3 0.3 0.3]*sz;

V = zeros(sz,sz,sz);
[ind1,ind2,ind3] = ndgrid(1:sz);
ind1 = ind1 - (sz+1)/2;
ind2 = ind2 - (sz+1)/2;
ind3 = ind3 - (sz+1)/2;

indA = [ind1(:) ind2(:) ind3(:)] * A;
inda = reshape(indA(:,1),size(ind1));
indb = reshape(indA(:,2),size(ind1));
indc = reshape(indA(:,3),size(ind1));

for i = 1:length(ctr1)
    r1 = sqrt((indb-ctr1(i)).^2 + (indc-ctr2(i)).^2);
    V10 = exp(-(r1/d1(i)).^2);
    V1 = exp(-((indb-ctr1(i))/d1(i)).^2 - ((indc-ctr2(i))/d2(i)).^2);
    V = V + V1;
end

[C,ah,bh] = get_volume_angle(V,[],'method','svd');

ahn = ah/norm(ah);
orthoplot(ind1,ind2,ind3,V, [0;0],[0;0],[0;0],[a(1);ahn(1)],[a(2);ahn(2)],[a(3);ahn(3)]);

fprintf('Actual a vector: [');
fprintf('%f ',a);
fprintf(']\n');

fprintf('Estimated a vector: [');
fprintf('%f ',ahn);
fprintf(']\n');

fprintf('Dot product: %f\n', dot(a,ahn));

fprintf('Actual b vector: [');
fprintf('%f ',b);
fprintf(']\n');

fprintf('Estimated b vector: [');
fprintf('%f ',bh);
fprintf(']\n');

fprintf('Dot product: %f\n', dot(b,bh));




