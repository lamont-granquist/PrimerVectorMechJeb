using System;
using System.Collections.Generic;
using System.Threading;
using UnityEngine;

namespace MuMech {
    public class PontryaginNode : PontryaginBase {
        public PontryaginNode(double mu, Vector3d r0, Vector3d v0, Vector3d pv0, Vector3d pr0, double dV, double bt) : base(mu, r0, v0, pv0, pr0, dV)
        {
            tc1 = bt / 2.0;
            tgo = bt;
            tc2 = 0;
            tc1_bar = tc1 / t_scale;
            tgo_bar = tgo / t_scale;
            tc2_bar = tc2 / t_scale;

            fixed_final_time = true;
        }

        private Orbit target;
        private double tc1, tc1_bar;
        private double tc2, tc2_bar;
        private double Tf;
        private Vector3d rT;
        private Vector3d vT;

        private bool initialized = false;

        public void intercept(Orbit target)
        {
            this.target = target;
            if (!initialized)
            {
                bcfun = intercept; // terminal5constraint;
                initialized = true;
            }
        }

        private void intercept(double[] yT, double[] z)
        {
            z[0] = yT[0] - rT[0];
            z[1] = yT[1] - rT[1];
            z[2] = yT[2] - rT[2];
            z[3] = yT[3] - vT[0];
            z[4] = yT[4] - vT[1];
            z[5] = yT[5] - vT[2];
        }

        private void terminal5constraint(double[] yT, double[] z)
        {
            Vector3d rTp = rT;
            Vector3d vTp = vT;
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d hT = Vector3d.Cross(rTp, vTp);

            if (Math.Abs(hT[1]) <= 1e-6)
            {
                // FIXME: this needs to go to temporary variables or needs to unfuck itself afterwards so that rf/vf can be updated as they're external
                rTp = rTp.Reorder(231);
                vTp = vTp.Reorder(231);
                rf = rf.Reorder(231);
                vf = vf.Reorder(231);
                prf = prf.Reorder(231);
                pvf = pvf.Reorder(231);
                hT = Vector3d.Cross(rTp, vTp);
            }

            Vector3d hf = Vector3d.Cross(rf, vf);
            Vector3d eT = - ( rTp.normalized + Vector3d.Cross(hT, vTp) );
            Vector3d ef = - ( rf.normalized + Vector3d.Cross(hf, vf) );
            Vector3d hmiss = hf - hT;
            Vector3d emiss = ef - eT;
            double trans = Vector3d.Dot(prf, vf) - Vector3d.Dot(pvf, rf) / ( rf.magnitude * rf.magnitude * rf.magnitude );

            z[0] = hmiss[0];
            z[1] = hmiss[1];
            z[2] = hmiss[2];
            z[3] = emiss[0];
            z[4] = emiss[2];
            z[5] = trans;
        }

        public override void Bootstrap(double t0)
        {
            // set the final time guess
            Tf = t0 + tgo * 3.0 / 2.0;
            Vector3d rf, vf;

            target.GetOrbitalStateVectorsAtUT(Tf, out rf, out vf);
            QuaternionD rot = Quaternion.Inverse(Planetarium.fetch.rotation);
            rT = rot * (rf.xzy / r_scale);
            vT = rot * (vf.xzy / v_scale);

            // build arcs off of ksp stages, with coasts
            List<Arc> arcs = new List<Arc>();

            arcs.Add(new Arc(new Stage(this, m0: stages[0].m0, isp: 0, thrust: 0, ksp_stage: stages[0].ksp_stage)));

            for(int i = 0; i < stages.Count; i++)
            {
                arcs.Add(new Arc(stages[i]));
            }

            arcs.Add(new Arc(new Stage(this, m0: stages[stages.Count-1].m0, isp: 0, thrust: 0, ksp_stage: stages[stages.Count-1].ksp_stage), done: true));
            // arcs.Add(new Arc(new Stage(this, m0: -1, isp: 0, thrust: 0, ksp_stage: stages[stages.Count-1].ksp_stage), done: true));

            arcs[arcs.Count-1].use_fixed_time = true;
            arcs[arcs.Count-1].fixed_time = Tf;
            arcs[arcs.Count-1].fixed_tbar = ( Tf - t0 ) / t_scale;

            // allocate y0
            y0 = new double[arcIndex(arcs, arcs.Count)];

            // initialize overall burn time (FIXME: duplicates code with UpdateY0)
            y0[0] = tgo_bar;
            y0[1] = 0;

            UpdateY0(arcs);

            // add guesses for coast burntimes
            y0[arcIndex(arcs, 0, parameters: true)] = tc1_bar;
            y0[arcIndex(arcs, 0, parameters: true)+1] = 0;

            // seed continuity initial conditions
            yf = new double[arcs.Count*13];
            multipleIntegrate(y0, yf, arcs, initialize: true);

            double[] y0_saved = new double[y0.Length];
            Array.Copy(y0, y0_saved, y0.Length);

            //for(int j = 0; j < y0.Length; j++)
            //    Debug.Log("bootstrap - y0[" + j + "] = " + y0[j]);
            //Debug.Log("running optimizer");

            bool success = runOptimizer(arcs);

            if ( !success )
            {
                //for(int k = 0; k < y0.Length; k++)
                //    Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                //Debug.Log("optimizer failed12");
                y0 = null;
                return;
            }

            Solution new_sol = new Solution(t_scale, v_scale, r_scale, t0);
            multipleIntegrate(y0, new_sol, arcs, 10);

            //Debug.Log("optimizer done");
            if ( new_sol.tgo(new_sol.t0, 0) < 0 )
            {
                // coast is less than zero, delete it and reconverge
                RemoveArc(arcs, 0, new_sol);
                UpdateY0(arcs); // reset to current starting position
                multipleIntegrate(y0, yf, arcs, initialize: true);
                //Debug.Log("running optimizer4");
                bcfun = terminal5constraint;

                if ( !runOptimizer(arcs) )
                {
                    //for(int k = 0; k < y0.Length; k++)
                    //    Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                    //Debug.Log("optimizer failed14");
                    y0 = null;
                    return;
                }

                //Debug.Log("optimizer done");
                new_sol = new Solution(t_scale, v_scale, r_scale, t0);
                multipleIntegrate(y0, new_sol, arcs, 10);
            }

            //for(int k = 0; k < y0.Length; k++)
                //Debug.Log("new y0[" + k + "] = " + y0[k]);

            this.solution = new_sol;
            //Debug.Log("done with bootstrap");
        }
    }
}
