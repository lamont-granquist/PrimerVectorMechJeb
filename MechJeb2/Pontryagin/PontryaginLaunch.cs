using System;
using System.Collections.Generic;
using System.Threading;
using UnityEngine;

namespace MuMech {
    public class PontryaginLaunch : PontryaginBase {
        public PontryaginLaunch(MechJebCore core, double mu, Vector3d r0, Vector3d v0, Vector3d pv0, Vector3d pr0, double dVguess = 0, double dTguess = 0) : base(core, mu, r0, v0, pv0, pr0)
        {
            this.dVguess = dVguess;
            this.dTguess = dTguess;
        }

        double dVguess;
        double dTguess;

        public bool omitCoast;

        double rTm;
        double vTm;
        double gammaT;
        double inc;
        double smaT;
        int numStages;
        Vector3d hT;

        // 5-constraint PEG with fixed LAN
        public void flightangle5constraint(double rTm, double vTm, double gammaT, Vector3d hT)
        {
            QuaternionD rot = Quaternion.Inverse(Planetarium.fetch.rotation);
            this.rTm = rTm / r_scale;
            this.vTm = vTm / v_scale;
            this.gammaT = gammaT;
            this.hT = rot * hT / r_scale / v_scale;
            bcfun = flightangle5constraint;
        }

        private void flightangle5constraint(double[] yT, double[] z)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);
            Vector3d hf = Vector3d.Cross(rf, vf);

            // rot takes hT and rotates it to the south pole ([ 0, -1, 0 ] in xzy coordinates in KSP)
            Quaternion rot = Quaternion.FromToRotation(hT, Vector3d.down);
            // now rotate our projected final vector by the same amount
            Vector3d hfp = rot * hf.normalized;

            /* 5 constraints */
            z[0] = ( rf.magnitude * rf.magnitude - rTm * rTm ) / 2.0;
            z[1] = ( vf.magnitude * vf.magnitude - vTm * vTm ) / 2.0;
            z[2] = Vector3d.Dot(rf, vf) - rf.magnitude * vf.magnitude * Math.Sin(gammaT);

            // stereographic projection onto the Riemann sphere ( hfp[1] is z-axis in KSP )
            z[3] = hfp[0] / ( 1.0 - hfp[1] );
            z[4] = hfp[2] / ( 1.0 - hfp[1] );

            // for an almost precisely 180 degree flip of the orbital angular momentum vector remap the infinity to something the
            // least squares algorithm can digest.
            if ( !z[3].IsFinite() || !z[4].IsFinite() || Math.Abs(z[3]) > Math.Sqrt(Double.MaxValue) || Math.Abs(z[4]) > Math.Sqrt(Double.MaxValue) )
            {
                z[3] = Math.Sqrt(Double.MaxValue);
                z[4] = Math.Sqrt(Double.MaxValue);
            }

            /* transversality - free argp */
            z[5] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), hT);
        }

        // 4-constraint PEG with free LAN
        public void flightangle4constraint(double rTm, double vTm, double gammaT, double inc)
        {
            //Debug.Log("call stack: + " + Environment.StackTrace);
            this.rTm = rTm / r_scale;
            this.vTm = vTm / v_scale;
            //Debug.Log("4constraint vTm = " + vTm + " v_scale = " + v_scale + " vTm_bar = " + this.vTm );
            //Debug.Log("4constraint rTm = " + rTm + " r_scale = " + r_scale + " rTm_bar = " + this.rTm );
            this.gammaT = gammaT;
            this.inc = inc;
            bcfun = flightangle4constraint;
        }

        private void flightangle4constraint(double[] yT, double[] z)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d n = new Vector3d(0, -1, 0);  /* angular momentum vectors point south in KSP and we're in xzy coords */
            Vector3d rn = Vector3d.Cross(rf, n);
            Vector3d vn = Vector3d.Cross(vf, n);
            Vector3d hf = Vector3d.Cross(rf, vf);

            z[0] = ( rf.magnitude * rf.magnitude - rTm * rTm ) / 2.0;
            z[1] = ( vf.magnitude * vf.magnitude - vTm * vTm ) / 2.0;
            z[2] = Vector3d.Dot(n, hf) - hf.magnitude * Math.Cos(inc);
            z[3] = Vector3d.Dot(rf, vf) - rf.magnitude * vf.magnitude * Math.Sin(gammaT);
            z[4] = rTm * rTm * ( Vector3d.Dot(vf, prf) - vTm * Math.Sin(gammaT) / rTm * Vector3d.Dot(rf, prf) ) -
                vTm * vTm * ( Vector3d.Dot(rf, pvf) - rTm * Math.Sin(gammaT) / vTm * Vector3d.Dot(vf, pvf) );
            z[5] = Vector3d.Dot(hf, prf) * Vector3d.Dot(hf, rn) + Vector3d.Dot(hf, pvf) * Vector3d.Dot(hf, vn);
        }

        // suborbital 3-constraint with fixed burntime and maximizes velocity
        public void suborbital3constraint(double rT, double gammaT, int numStages, double inc)
        {
            this.rTm = rTm / r_scale;
            this.gammaT = gammaT;
            this.numStages = numStages;
            this.inc = inc;
            bcfun = suborbital3constraint;
        }

        private void suborbital3constraint(double[] yT, double[] z)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d n = new Vector3d(0, -1, 0);  /* angular momentum vectors point south in KSP and we're in xzy coords */
            Vector3d rn = Vector3d.Cross(rf, n);
            Vector3d vn = Vector3d.Cross(vf, n);
            Vector3d hf = Vector3d.Cross(rf, vf);

            z[0] = ( rf.magnitude * rf.magnitude - rTm * rTm ) / 2.0;
            z[1] = Vector3d.Dot(n, Vector3d.Cross(rf, vf)) - Vector3d.Cross(rf, vf).magnitude * Math.Cos(inc);
            z[2] = Vector3d.Dot(rf, vf) - rf.magnitude * vf.magnitude * Math.Sin(gammaT);

            z[3] = Vector3d.Dot(vf, prf) * rf.sqrMagnitude - Vector3d.Dot(rf, pvf) * vf.sqrMagnitude
                + Vector3d.Dot(rf, vf) * (vf.sqrMagnitude - Vector3d.Dot(rf, prf));
            z[4] = Vector3d.Dot(vf, pvf) - vf.sqrMagnitude;
            z[5] = Vector3d.Dot(hf, prf) * Vector3d.Dot(hf, Vector3d.Cross(rf, n))
                + Vector3d.Dot(hf, pvf) * Vector3d.Dot(hf, Vector3d.Cross(vf, n));
        }

        public override void Bootstrap(double t0)
        {
            // build arcs off of ksp stages, with coasts
            List<Arc> arcs = new List<Arc>();
            Debug.Log("STAGES BEFORE BOOTSTRAP:");
            for(int i = 0; i < stages.Count; i++)
            {
                Debug.Log(stages[i]);
                arcs.Add(new Arc(this, stages[i]));
            }

            arcs[arcs.Count-1].infinite = true;

            // allocate y0
            y0 = new double[arcIndex(arcs, arcs.Count)];

            // halfway wild guess at burntime to orbit
            if (dVguess > 0)
                tgo = stages[0].v_e * stages[0].startMass / stages[0].startThrust * ( 1 - Math.Exp(-dVguess/stages[0].v_e) );
            else if (dTguess > 0)
                tgo = dTguess;
            else
                tgo = 240;

            tgo_bar = tgo / t_scale;

            // initialize overall burn time
            y0[0] = tgo_bar;
            y0[1] = 0;

            UpdateY0(arcs);

            // seed continuity initial conditions
            yf = new double[arcs.Count*13];
            multipleIntegrate(y0, yf, arcs, initialize: true);

            //for(int j = 0; j < y0.Length; j++)
            //    Debug.Log("bootstrap - y0[" + j + "] = " + y0[j]);

            //Debug.Log("running optimizer");

            if ( !runOptimizer(arcs) )
            {
                //for(int k = 0; k < y0.Length; k++)
                    //Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                //Debug.Log("optimizer failed5");
                y0 = null;
                return;
            }

            //Debug.Log("optimizer done");

            Solution new_sol = new Solution(t_scale, v_scale, r_scale, t0);
            multipleIntegrate(y0, new_sol, arcs, 10);

            for(int i = 0; i < arcs.Count; i++)
            {
                if ( new_sol.tgo(new_sol.t0, i) < 1 )
                {
                    /* burn is less than one second, delete it */
                    RemoveArc(arcs, i, new_sol);
                }
            }

            //for(int k = 0; k < y0.Length; k++)
            //   Debug.Log("y0[" + k + "] = " + y0[k]);

            //Debug.Log("arcs.Count = " + arcs.Count);
            //Debug.Log("segments.Count = " + new_sol.segments.Count);

            bool insertedCoast = false;

            if (arcs.Count > 1 && !omitCoast)
            {
                InsertCoast(arcs, arcs.Count-1, new_sol);
                insertedCoast = true;
            }
            arcs[arcs.Count-1].infinite = false;

            //Debug.Log("running optimizer3");

            if ( !runOptimizer(arcs) )
            {
                //for(int k = 0; k < y0.Length; k++)
                    //Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                //Debug.Log("optimizer failed6");
                y0 = null;
                return;
            }

            //Debug.Log("optimizer done");

            new_sol = new Solution(t_scale, v_scale, r_scale, t0);
            multipleIntegrate(y0, new_sol, arcs, 10);

            //for(int k = 0; k < y0.Length; k++)
            //    Debug.Log("y0[" + k + "] = " + y0[k]);

            if (insertedCoast)
            {

                /*
                if ( new_sol.tgo(new_sol.t0, arcs.Count-2) < 1 )
                {
                    // coast is less than one second, try extending it again without infinite
                    // (FIXME: this is buggy somehow)
                    RetryCoast(arcs, arcs.Count-2, new_sol);
                    //Debug.Log("running optimizer4");

                    if ( !runOptimizer(arcs) )
                    {
                        for(int k = 0; k < y0.Length; k++)
                            Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                        Debug.Log("optimizer failed8");
                        y0 = null;
                        return;
                    }

                    //Debug.Log("optimizer done");
                    new_sol = new Solution(t_scale, v_scale, r_scale, t0);
                    multipleIntegrate(y0, new_sol, arcs, 10);
                }
                */

                if ( new_sol.tgo(new_sol.t0, arcs.Count-2) < 1 )
                {
                    // coast is less than one second, delete it and reconverge
                    RemoveArc(arcs, arcs.Count-2, new_sol);
                    //Debug.Log("running optimizer4");

                    if ( !runOptimizer(arcs) )
                    {
                        //for(int k = 0; k < y0.Length; k++)
                        //    Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                        //Debug.Log("optimizer failed10");
                        y0 = null;
                        return;
                    }

                    //Debug.Log("optimizer done");
                    new_sol = new Solution(t_scale, v_scale, r_scale, t0);
                    multipleIntegrate(y0, new_sol, arcs, 10);
                }
            }

            //for(int k = 0; k < y0.Length; k++)
                //Debug.Log("new y0[" + k + "] = " + y0[k]);

            this.solution = new_sol;
            //Debug.Log("done with bootstrap");

            yf = new double[arcs.Count*13];
            multipleIntegrate(y0, yf, arcs);

            //for(int k = 0; k < yf.Length; k++)
                //Debug.Log("new yf[" + k + "] = " + yf[k]);

            //Debug.Log("optimizer hT = " + hT.magnitude * r_scale * v_scale);
            //Debug.Log("r_scale = " + r_scale);
            //Debug.Log("v_scale = " + v_scale);
            //Debug.Log("rTm = " + rTm * r_scale + " " + rTm);
            //Debug.Log("vTm = " + vTm * v_scale + " " + vTm);
            Debug.Log("STAGES AFTER BOOTSTRAP:");
            for(int i = 0; i < stages.Count; i++)
            {
                Debug.Log(stages[i]);
            }
            Debug.Log("ARCS AFTER BOOTSTRAP:");
            for(int i = 0; i < arcs.Count; i++)
            {
                Debug.Log(arcs[i]);
            }
        }

        public override void Optimize(double t0)
        {
            base.Optimize(t0);
            //Debug.Log("r_scale = " + r_scale);
            //Debug.Log("v_scale = " + v_scale);
            //Debug.Log("rTm = " + rTm * r_scale + " " + rTm);
            //Debug.Log("vTm = " + vTm * v_scale + " " + vTm);
        }
    }
}
