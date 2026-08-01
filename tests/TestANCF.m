classdef TestANCF < matlab.unittest.TestCase
    methods (Test)
        function straightElementHasNegligibleInternalForce(testCase)
            p.Le = 0.25; p.E = 2.1e9; p.rho = 7800;
            p.A = 2e-5; p.I = 1.67e-11; p.g = 9.81;
            q = [0;0;1;0;p.Le;0;1;0];
            [M,fg,fi] = ancfBeam2D(q,p);
            testCase.verifyLessThan(norm(fi),1e-8);
            testCase.verifyEqual(M,M.','AbsTol',1e-12);
            testCase.verifyGreaterThan(min(eig((M+M.')/2)),0);
            testCase.verifyLessThan(fg(2)+fg(6),0);
        end

        function interpolationReproducesEndCoordinates(testCase)
            L = 0.4;
            q = [0.1;-0.2;0.9;0.1;0.45;0.15;0.8;0.2];
            S0 = ancfShape2D(0,L);
            S1 = ancfShape2D(1,L);
            testCase.verifyEqual(S0*q,q(1:2),'AbsTol',1e-14);
            testCase.verifyEqual(S1*q,q(5:6),'AbsTol',1e-14);
        end

        function flexiblePendulumSmokeTest(testCase)
            out = run_ancf_flexible_pendulum(false);
            testCase.verifyTrue(all(isfinite(out.q),'all'));
            testCase.verifyEqual(out.t(end),0.35,'AbsTol',1e-10);
            testCase.verifyLessThan(min(out.tip(:,2)),-1e-4);
            testCase.verifyLessThan(abs(out.q(end,1)),1e-14);
            testCase.verifyLessThan(abs(out.q(end,2)),1e-14);
        end
    end
end
