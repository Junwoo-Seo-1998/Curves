using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Reflection.Emit;
using System.Security;
using System.Security.Cryptography;
using System.Text;
using System.Windows.Forms;
using Label = System.Windows.Forms.Label;

namespace mat_290_framework
{
    public partial class MAT290 : Form
    {
        private List<List<int>> PascalValues;
        private List<List<Point2D>> NewtonFormList;
        private List<double> splineXCoeff;
        private List<double> splineYCoeff;

        private float padding_left = 0;
        private float padding_right = 0;

        private float padding_top = 0;
        private float padding_bottom = 0;

        private float graph_y_min_max = 5;
        private float graph_x_min_max = 10;

        private float y_min_max = 3.0f;
        const float rad = 25.0f;
        private List<Label> labels=new List<Label>();
        private List<TrackBar> trackBars = new List<TrackBar>();
        private void MakeLabel(string text, Point2D pos)
        {
            Label templab = new Label();
            templab.Text = text;
            templab.Location = NonNegativeGraphPointToWindowPoint(pos).P();
            //templab.AutoSize = true;
            templab.Font= new Font("Arial", 10);
            templab.BackColor = System.Drawing.Color.Transparent;
            // Adding this control to the form
            Controls.Add(templab);
            labels.Add(templab);
        }

        private void RemoveLabels()
        {
            foreach (var label in labels)
            {
                Controls.Remove(label);
            }
        }
        void OnTrackBarValueChanged(object sender, EventArgs e)
        {
            // get trackbar, which generated event
            var trackBar = (TrackBar)sender;
            float value = trackBar.Value / 100.0f;
            tVal_ = value;
            Refresh();
        }
        private void MakeTvalSlider()
        {
            TrackBar trackBar = new TrackBar();
            trackBar.ValueChanged += OnTrackBarValueChanged;
            trackBar.Minimum = 0;
            trackBar.Maximum = 100;
            trackBar.Value = 50;
            trackBar.Width = 200;
            int w = Bounds.Width;
            int h = Bounds.Height;
            trackBar.Location = new Point(w - 250, h - 100);
            Controls.Add(trackBar);
            trackBars.Add(trackBar);
        }

        private void RemoveSlider()
        {
            foreach (var tbar in trackBars)
            {
                Controls.Remove(tbar);
            }
        }

        public MAT290()
        {
            InitializeComponent();
            padding_left = Bounds.Width * 0.1f;
            padding_right = Bounds.Width * 0.1f;
            padding_top = Bounds.Height * 0.1f;
            padding_bottom = Bounds.Height * 0.1f;
            pts_ = new List<Point2D>();
            tVal_ = 0.5F;
            degree_ = 0;
            knot_ = new List<float>();
            EdPtCont_ = true;
            rnd_ = new Random();
            PascalValues=new List<List<int>>();
            NewtonFormList = new List<List<Point2D>>();
            splineXCoeff = new List<double>();
            splineYCoeff = new List<double>();

            for (int i = 0; i < 21; ++i)
            {
                PascalValues.Add(new List<int>());
            }

            PascalValues[0].Add(1);
            for (int row = 1; row < 21; ++row)
            {
                PascalValues[row].Add(1);
                for (int i = 1; i < row; ++i)
                {
                    PascalValues[row].Add(PascalValues[row - 1][i - 1] + PascalValues[row - 1][i]);
                }
                PascalValues[row].Add(1);
            }

            shellLinePen.DashPattern = dashValues;
            Menu_DeBoor_Click(null, null);
        }

        // Point class for general math use
        protected class Point2D : System.Object
        {
            public double x;
            public double y;

            public Point2D(float _x, float _y)
            {
                x = _x;
                y = _y;
            }

            public Point2D(double _x, double _y)
            {
                x = _x;
                y = _y;
            }

            public Point2D(Point2D rhs)
            {
                x = rhs.x;
                y = rhs.y;
            }

            // adds two points together; used for barycentric combos
            public static Point2D operator +(Point2D lhs, Point2D rhs)
            {
                return new Point2D(lhs.x + rhs.x, lhs.y + rhs.y);
            }

            public static Point2D operator -(Point2D lhs, Point2D rhs)
            {
                return new Point2D(lhs.x - rhs.x, lhs.y - rhs.y);
            }

            public static Point2D operator /(Point2D lhs, float val)
            {
                return new Point2D(lhs.x / val, lhs.y / val);
            }

            public static Point2D operator /(Point2D lhs, double val)
            {
                return new Point2D(lhs.x / val, lhs.y / val);
            }

            // gets a distance between two points. not actual distance; used for picking
            public static float operator %(Point2D lhs, Point2D rhs)
            {
                float dx = (float)(lhs.x - rhs.x);
                float dy = (float)(lhs.y - rhs.y);

                return (dx * dx + dy * dy);
            }

            // scalar multiplication of points; for barycentric combos
            public static Point2D operator *(float t, Point2D rhs)
            {
                return new Point2D(rhs.x * t, rhs.y * t);
            }

            // scalar multiplication of points; for barycentric combos
            public static Point2D operator *(Point2D rhs, float t)
            {
                return new Point2D(rhs.x * t, rhs.y * t);
            }

            public static Point2D operator *(Point2D rhs, double t)
            {
                return new Point2D(rhs.x * t, rhs.y * t);
            }


            // returns the drawing subsytems' version of a point for drawing.
            public System.Drawing.Point P()
            {
                return new System.Drawing.Point((int)x, (int)y);
            }
        };

        List<Point2D> pts_; // the list of points used in internal algthms
        float tVal_; // t-value used for shell drawing
        int degree_; // degree of deboor subsplines
        List<float> knot_; // knot sequence for deboor
        bool EdPtCont_; // end point continuity flag for std knot seq contruction
        Random rnd_; // random number generator

        // pickpt returns an index of the closest point to the passed in point
        //  -- usually a mouse position
        private int PickPt(Point2D m)
        {
            float closest = m % pts_[0];
            int closestIndex = 0;

            for (int i = 1; i < pts_.Count; ++i)
            {
                float dist = m % pts_[i];
                if (dist < closest)
                {
                    closest = dist;
                    closestIndex = i;
                }
            }

            return closestIndex;
        }

        private void Menu_Clear_Click(object sender, EventArgs e)
        {
            pts_.Clear();
            Refresh();
        }

        private void Menu_Exit_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        private void MAT290_MouseMove(object sender, MouseEventArgs e)
        {
            if (Menu_DeCast.Checked || Menu_Bern.Checked)
            {
                if (pts_.Count != 0 && e.Button == MouseButtons.Left)
                {
                    // grab the closest point and snap it to the mouse
                    int index = PickPt(new Point2D(e.X, e.Y));

                    pts_[index].y = e.Y;

                    Refresh();
                }
                return;
            }

            // if the right mouse button is being pressed
            if (pts_.Count != 0 && e.Button == MouseButtons.Right)
            {
                // grab the closest point and snap it to the mouse
                int index = PickPt(new Point2D(e.X, e.Y));

                pts_[index].x = e.X;
                pts_[index].y = e.Y;

                Refresh();
            }
        }

        private void MAT290_MouseDown(object sender, MouseEventArgs e)
        {
            int max_points = 20;

            if (Menu_DeBoor.Checked)
            {
                max_points = 40;
            }

            if (pts_.Count == max_points)
                    return;

            if (Menu_DeCast.Checked || Menu_Bern.Checked)
            {
                return;
            }

            // if the left mouse button was clicked
            if (e.Button == MouseButtons.Left)
            {
                
                // add a new point to the controlPoints
                pts_.Add(new Point2D(e.X, e.Y));

                if (Menu_DeBoor.Checked)
                {
                    ResetKnotSeq();
                    UpdateKnotSeq();
                }
                else
                {
                    degree_ = pts_.Count - 1;
                }

                Refresh();
            }

            // if there are points and the middle mouse button was pressed
            if (pts_.Count != 0 && e.Button == MouseButtons.Middle)
            {
                // then delete the closest point
                int index = PickPt(new Point2D(e.X, e.Y));

                pts_.RemoveAt(index);

                if (Menu_DeBoor.Checked)
                {
                    ResetKnotSeq();
                    UpdateKnotSeq();
                }

                Refresh();
            }
        }

        private void MAT290_MouseWheel(object sender, MouseEventArgs e)
        {
            // if the mouse wheel has moved
            if (e.Delta != 0)
            {
                // change the t-value for shell
                tVal_ += e.Delta / 120 * .02f;

                // handle edge cases
                tVal_ = (tVal_ < 0) ? 0 : tVal_;
                tVal_ = (tVal_ > 1) ? 1 : tVal_;

                Refresh();
            }
        }

        private void ResetPoints()
        {
            int num_of_points = degree_ + 1;
            float space = 1.0f / degree_;
            pts_.Clear();


            for (int i = 0; i < num_of_points; ++i)
            {
                Point2D p = NonNegativeGraphPointToWindowPoint(new Point2D(i * space, 1.0f));
                pts_.Add(p);
            }
        }

        private void NUD_degree_ValueChanged(object sender, EventArgs e)
        {
            if (pts_.Count == 0)
                return;

            degree_ = (int)NUD_degree.Value;
            if (Menu_DeCast.Checked || Menu_Bern.Checked)
            {
                //clamp
                degree_ = (degree_ < 1) ? 1 : degree_;
                degree_ = (degree_ > 20) ? 20 : degree_;
                NUD_degree.Value = degree_;
                ResetPoints();
            }
            else
            {
                ResetKnotSeq();
                UpdateKnotSeq();

                NUD_degree.Value = degree_;
            }
            Refresh();
        }

        private void CB_cont_CheckedChanged(object sender, EventArgs e)
        {
            EdPtCont_ = CB_cont.Checked;

            ResetKnotSeq();
            UpdateKnotSeq();

            Refresh();
        }

        private void Txt_knot_KeyPress(object sender, KeyPressEventArgs e)
        {
            if (e.KeyChar == '\r' || e.KeyChar == '\n')
            {
                // update knot seq
                string[] splits = Txt_knot.Text.ToString().Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

                if (splits.Length > pts_.Count + degree_ + 1)
                    return;

                knot_.Clear();
                foreach (string split in splits)
                {
                    knot_.Add(Convert.ToSingle(split));
                }

                for (int i = knot_.Count; i < (pts_.Count + degree_ + 1); ++i)
                    knot_.Add((float)(i - degree_));

                UpdateKnotSeq();
            }

            Refresh();
        }

        private void ResetMenus()
        {
            Menu_BezierCurves_DeCast.Checked = false;
            Menu_BezierCurves_Bern.Checked = false;

            Menu_Bern.Checked = false;
            Menu_DeCast.Checked = false;

            Menu_Midpoint.Checked = false;

            Menu_Inter_Poly.Checked = false;
            Menu_Inter_Splines.Checked = false;

            Menu_DeBoor.Checked = false;

            Menu_Polyline.Enabled = true;
            Menu_Points.Enabled = true;
            Menu_Shell.Enabled = false;
            Menu_Shell.Checked=false;

            pts_.Clear();
            RemoveSlider();
        }
        private void Menu_Polyline_Click(object sender, EventArgs e)
        {
            Refresh();
        }

        private void Menu_Points_Click(object sender, EventArgs e)
        {
            Refresh();
        }

        private void Menu_Shell_Click(object sender, EventArgs e)
        {
            Refresh();
        }

        private void Menu_DeCast_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_DeCast.Checked = true;

            ToggleDeBoorHUD(false);
            Lbl_degree.Visible = true;
            NUD_degree.Visible = true;

            degree_ = 3;
            ResetPoints();
            NUD_degree.Value = degree_;
            
            Refresh();
        }

        private void Menu_Bern_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_Bern.Checked = true;

            ToggleDeBoorHUD(false);
            Lbl_degree.Visible = true;
            NUD_degree.Visible = true;

            degree_ = 3;
            ResetPoints();
            NUD_degree.Value = degree_;

            Refresh();
        }

        private void Menu_Bezier_DeCast_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_BezierCurves_DeCast.Checked = true;
            Menu_Shell.Enabled = true;
            Menu_Shell.Checked = true;
            ToggleDeBoorHUD(false);
            Lbl_degree.Visible = false;
            NUD_degree.Visible = false;
            MakeTvalSlider();

            Refresh();
        }

        private void Menu_BezierCurves_Bern_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_BezierCurves_Bern.Checked = true;

            ToggleDeBoorHUD(false);
            Lbl_degree.Visible = false;
            NUD_degree.Visible = false;


            Refresh();
        }

        private void Menu_Midpoint_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_Midpoint.Checked = true;
            ToggleDeBoorHUD(false);
            Refresh();
        }

        private void Menu_Inter_Poly_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_Inter_Poly.Checked = true;

            ToggleDeBoorHUD(false);

            Refresh();
        }

        private void Menu_Inter_Splines_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_Inter_Splines.Checked = true;

            ToggleDeBoorHUD(false);
            RemoveLabels();

            Refresh();
        }

        private void Menu_DeBoor_Click(object sender, EventArgs e)
        {
            ResetMenus();
            Menu_DeBoor.Checked = true;

            ToggleDeBoorHUD(true);
            RemoveLabels();

            Refresh();
        }

        private void DegreeClamp()
        {
            // handle edge cases
            degree_ = (degree_ > pts_.Count - 1) ? pts_.Count - 1 : degree_;
            degree_ = (degree_ < 1) ? 1 : degree_;
            //clamp max
            degree_ = (degree_ > 20) ? 20 : degree_;
        }

        private void ResetKnotSeq( )
        {
            DegreeClamp();
            knot_.Clear();

            if (EdPtCont_)
            {
                for (int i = 0; i < degree_; ++i)
                    knot_.Add(0.0f);
                for (int i = 0; i <= (pts_.Count - degree_); ++i)
                    knot_.Add((float)i);
                for (int i = 0; i < degree_; ++i)
                    knot_.Add((float)(pts_.Count - degree_));
            }
            else
            {
                for (int i = -degree_; i <= (pts_.Count); ++i)
                    knot_.Add((float)i);
            }
        }

        private void UpdateKnotSeq()
        {
            Txt_knot.Clear();
            foreach (float knot in knot_)
            {
                Txt_knot.Text += knot.ToString() + " ";
            }
        }

        private void ToggleDeBoorHUD( bool on )
        {
            // set up basic knot sequence
            if(on)
            {
                ResetKnotSeq();
                UpdateKnotSeq();
            }

            CB_cont.Visible = on;

            Lbl_knot.Visible = on;
            Txt_knot.Visible = on;

            Lbl_degree.Visible = on;
            NUD_degree.Visible = on;

            Lbl_N.Visible = on;
            NUD_N.Visible = on;
        }

        private void MAT290_Paint(object sender, PaintEventArgs e)
        {
            // pass the graphics object to the DrawScreen subroutine for processing
            DrawScreen(e.Graphics);
        }

        private Point2D WindowPointToGraphPoint(Point2D point)
        {
            double w = Bounds.Width - (padding_left + padding_right);
            double h = Bounds.Height - (padding_top + padding_bottom);

            double A = w / (2.0 * graph_x_min_max);
            double B = -h / (2.0 * graph_y_min_max);

            double tX = w / 2.0 + padding_left;
            double tY = h / 2.0 + padding_top;

            tX = -tX / A;
            tY = -tY / B;
            
            A = 1.0 / A;
            B = 1.0 / B;

            return new Point2D(A * point.x + tX, B * point.y + tY);
        }

        private Point2D GraphPointToWindowPoint(Point2D point)
        {
            double w = Bounds.Width - (padding_left + padding_right);
            double h = Bounds.Height - (padding_top + padding_bottom);

            double A = w / (2.0 * graph_x_min_max);
            double B = -h / (2.0 * graph_y_min_max);

            double tX = w / 2.0 + padding_left;
            double tY = h / 2.0 + padding_top;

            return new Point2D(A * point.x + tX, B * point.y + tY);
        }

        private Point2D WindowPointToNonNegativeGraphPoint(Point2D point)
        {
            double w = Bounds.Width;
            double h = Bounds.Height;

            double By = (h - padding_bottom + padding_top) / 2.0;
            double Ay = (padding_top - By) / (double)y_min_max;

            double x = point.x / (w - padding_right - padding_left) +
                       (-padding_left / (w - padding_right - padding_left));
            double y = (point.y / Ay) + (-By / Ay);
            return new Point2D((float)x, (float)y);
        }

        private Point2D NonNegativeGraphPointToWindowPoint(Point2D point)
        {
            double w = Bounds.Width;
            double h = Bounds.Height;

            double By = (h - padding_bottom + padding_top) / 2.0;
            double Ay = (padding_top - By) / (double)y_min_max;

            double x = (w - padding_right - padding_left) * point.x + padding_left;
            double y = Ay * point.y + By;

            return new Point2D((float)x, (float)y);
        }

        // pens used for drawing elements of the display
        System.Drawing.Pen polyPen = new Pen(Color.Gray, 1.0f);
        System.Drawing.Pen shellPen = new Pen(Color.Red, 0.5f);
        float[] dashValues = { 5, 2, 15, 4 };
        System.Drawing.Pen shellLinePen = new Pen(Color.Red, 0.5f);

        System.Drawing.Pen splinePen = new Pen(Color.Navy, 1.5f);

        System.Drawing.Pen xPen = new Pen(Color.Red, 0.5f);
        System.Drawing.Pen hlinePen = new Pen(Color.DimGray, 0.25f);
        System.Drawing.Pen yPen = new Pen(Color.DeepSkyBlue, 0.5f);
        Font bFont = new Font("Arial", 12);
        Font arial = new Font("Arial", 12);
        private int grid_step = 1;
        private void DrawGrid(System.Drawing.Graphics gfx)
        {
            
            var p1 = GraphPointToWindowPoint(new Point2D(0.0f, graph_y_min_max)).P();
            var p2 = GraphPointToWindowPoint(new Point2D(0.0f, -graph_y_min_max)).P();

            gfx.DrawLine(yPen, p1, p2);

            p1 = GraphPointToWindowPoint(new Point2D(-graph_x_min_max, 0)).P();
            p2 = GraphPointToWindowPoint(new Point2D(graph_x_min_max, 0)).P();
            gfx.DrawLine(xPen, p1, p2);

            p1 = GraphPointToWindowPoint(new Point2D(0, 0)).P();
            gfx.DrawString("0", bFont, Brushes.Red, p1.X, p1.Y);

            int steps = (int)graph_y_min_max / grid_step;
            for (int i = 1; i <= steps; i++)
            {
                float yval = i * grid_step;
                p1 = GraphPointToWindowPoint(new Point2D(-graph_x_min_max, yval)).P();
                p2 = GraphPointToWindowPoint(new Point2D(graph_x_min_max, yval)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                p1 = GraphPointToWindowPoint(new Point2D(0, yval)).P();
                gfx.DrawString(yval.ToString(), bFont, Brushes.Black, p1.X, p1.Y);

                yval = -yval;
                p1 = GraphPointToWindowPoint(new Point2D(-graph_x_min_max, yval)).P();
                p2 = GraphPointToWindowPoint(new Point2D(graph_x_min_max, yval)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                p1 = GraphPointToWindowPoint(new Point2D(0, yval)).P();
                gfx.DrawString(yval.ToString(), bFont, Brushes.Black, p1.X, p1.Y);
            }
            steps = (int)graph_x_min_max / grid_step;
            for (int i = 1; i <= steps; i++)
            {
                float xval = i * grid_step;
                p1 = GraphPointToWindowPoint(new Point2D(xval, graph_y_min_max)).P();
                p2 = GraphPointToWindowPoint(new Point2D(xval, -graph_y_min_max)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                p1 = GraphPointToWindowPoint(new Point2D(xval, 0)).P();
                gfx.DrawString(xval.ToString(), bFont, Brushes.Black, p1.X, p1.Y);

                xval = -xval;
                p1 = GraphPointToWindowPoint(new Point2D(xval, graph_y_min_max)).P();
                p2 = GraphPointToWindowPoint(new Point2D(xval, -graph_y_min_max)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                p1 = GraphPointToWindowPoint(new Point2D(xval, 0)).P();
                gfx.DrawString(xval.ToString(), bFont, Brushes.Black, p1.X, p1.Y);
            }
        }
        private void DrawScreen(System.Drawing.Graphics gfx)
        {
            if (Menu_DeCast.Checked)
            {
                gfx.DrawString("DeCasteljau", arial, Brushes.Black, 0, 30);
            }
            else if (Menu_BezierCurves_DeCast.Checked)
            {
                gfx.DrawString("DeCasteljau", arial, Brushes.Black, 0, 30);
            }
            else if (Menu_BezierCurves_Bern.Checked)
            {
                gfx.DrawString("Bernstein", arial, Brushes.Black, 0, 30);
            }
            else if (Menu_Midpoint.Checked)
            {
                gfx.DrawString("Midpoint", arial, Brushes.Black, 0, 30);
            }
            else if (Menu_Bern.Checked)
            {
                gfx.DrawString("Bernstein", arial, Brushes.Black, 0, 30);
            }
            else if (Menu_DeBoor.Checked)
            {
                DrawGrid(gfx);
                gfx.DrawString("DeBoor", arial, Brushes.Black, 0, 30);
            }

            // to prevent unecessary drawing
            if (pts_.Count == 0)
                return;

            if (Menu_Shell.Checked)
            {
                // draw the shell
                DrawShell(gfx, pts_, tVal_);
            }

            if (Menu_Polyline.Checked)
            {
                // draw the control poly
                for (int i = 1; i < pts_.Count; ++i)
                {
                    gfx.DrawLine(polyPen, pts_[i - 1].P(), pts_[i].P());
                }
            }

            if (Menu_Points.Checked)
            {
                // draw the control points
                foreach (Point2D pt in pts_)
                {
                    gfx.DrawEllipse(polyPen, (float)(pt.x - rad / 2.0), (float)(pt.y - rad / 2.0), rad, rad);
                    gfx.DrawEllipse(polyPen, (float)(pt.x - rad / 2.0), (float)(pt.y - rad / 2.0), rad, rad);
                }
            }

            // you can change these variables at will; i have just chosen there
            //  to be six sample points for every point placed on the screen
            float steps = pts_.Count * 6;
            float alpha = 1 / steps;

            ///////////////////////////////////////////////////////////////////////////////
            // Drawing code for algorithms goes in here                                  //
            ///////////////////////////////////////////////////////////////////////////////

            // DeCastlejau algorithm
            if (Menu_DeCast.Checked)
            {
                DrawDeCastlejau(gfx, alpha);
            }

            // Bernstein polynomial
            if (Menu_Bern.Checked)
            {
                DrawBernstein(gfx, alpha);
            }
            // BezierCurves using De Casteljau Algorithm
            if (Menu_BezierCurves_DeCast.Checked)
            {
                DrawBezier_DeCast(gfx, alpha);
            }

            if (Menu_BezierCurves_Bern.Checked)
            {
                DrawBezier_Bern(gfx, alpha);
            }

            // Midpoint algorithm
            if (Menu_Midpoint.Checked)
            {
                DrawMidpoint(gfx, pts_, 0);
            }

            // polygon interpolation
            if (Menu_Inter_Poly.Checked)
            {
                DrawPolyInterpolate(gfx, alpha);
            }

            // spline interpolation
            if (Menu_Inter_Splines.Checked)
            {
                DrawSplineInterpolate(gfx, alpha);
            }

            // deboor
            if (Menu_DeBoor.Checked && pts_.Count >= 2)
            {
                DrawDeBoor(gfx, alpha);
            }

            ///////////////////////////////////////////////////////////////////////////////
            // Drawing code end                                                          //
            ///////////////////////////////////////////////////////////////////////////////


            // Heads up Display drawing code

            if (Menu_BezierCurves_DeCast.Checked||Menu_Midpoint.Checked)
            {
                gfx.DrawString("t-value: " + tVal_.ToString("F"), arial, Brushes.Black, 500, 30);
            }

            //gfx.DrawString("t-step: " + alpha.ToString("F6"), arial, Brushes.Black, 600, 30);

            gfx.DrawString("points: "+pts_.Count.ToString(), arial, Brushes.Black, 750, 30);
        }

        private void DrawLineStrip(System.Drawing.Graphics gfx, System.Drawing.Pen pen, List<Point2D> pts)
        {
            if (pts.Count < 2)
                return;

            Point2D current_left = pts[0];
            Point2D current_right = pts[1];
            gfx.DrawLine(pen, current_left.P(), current_right.P());
            for (int i = 2; i < pts.Count; ++i)
            {
                current_left = current_right;
                current_right = pts[i];
                gfx.DrawLine(pen, current_left.P(), current_right.P());
            }
        }
        private void DrawShell(System.Drawing.Graphics gfx, List<Point2D> pts, float t)
        {
            List<Point2D> decast_vals = new List<Point2D>();
            List<Point2D> points_to_draw = new List<Point2D>();
            float one_min_t = 1.0f - t;
            for (int i = 0; i < pts.Count; ++i)
            {
                decast_vals.Add(pts[i]);
            }
            for (int deg = pts.Count - 1; deg >= 0; --deg)
            {
                points_to_draw.Clear();
                for (int i = 0; i < deg; ++i)
                {
                    decast_vals[i] = decast_vals[i] * one_min_t + decast_vals[i + 1] * t;
                    DrawCircle(gfx, shellPen, decast_vals[i]);
                    points_to_draw.Add(decast_vals[i]);
                }
                DrawLineStrip(gfx, shellLinePen, points_to_draw);
            }
        }

        private Point2D Gamma(int start, int end, float t)
        {
            return new Point2D(0, 0);
        }

        private Point2D DeCastlejau(float t)
        {
            List<float> decast_vals = new List<float>();
            float one_min_t = 1.0f - t;

            for (int i = 0; i < pts_.Count; ++i)
            {
                float y = (float)WindowPointToNonNegativeGraphPoint(pts_[i]).y;
                decast_vals.Add(y);
            }

            for (int deg = degree_; deg >= 0; --deg) 
            {
                for (int i = 0; i < deg; i++)
                {
                    decast_vals[i] = decast_vals[i] * one_min_t + decast_vals[i + 1] * t;
                }
            }

            var graphPointToWindow = NonNegativeGraphPointToWindowPoint(new Point2D(t, decast_vals[0]));
            return graphPointToWindow;
        }

        private Point2D Bernstein(float t)
        {
            List<float> one_min_t_pows=new List<float>();
            List<float> t_pows = new List<float>();
            one_min_t_pows.Add(1);
            t_pows.Add(1);
            for (int i = 1; i <= degree_; i++)
            {
                one_min_t_pows.Add(one_min_t_pows[i-1]*(1-t));
                t_pows.Add(t_pows[i-1] * t);
            }
            
            float sum = 0;
            for (int i = 0; i <= degree_; i++)
            {
                float c = (float)WindowPointToNonNegativeGraphPoint(new Point2D(0, pts_[i].y)).y;
                sum += c * PascalValues[degree_][i] * one_min_t_pows[degree_ - i] * t_pows[i];
            }
            var graphPointToWindow = NonNegativeGraphPointToWindowPoint(new Point2D(t, sum));
            return graphPointToWindow;
        }

        private Point2D Bezier_Bernstein(float t)
        {
            List<float> one_min_t_pows = new List<float>();
            List<float> t_pows = new List<float>();
            one_min_t_pows.Add(1);
            t_pows.Add(1);
            for (int i = 1; i <= degree_; i++)
            {
                one_min_t_pows.Add(one_min_t_pows[i - 1] * (1 - t));
                t_pows.Add(t_pows[i - 1] * t);
            }

            Point2D sum = new Point2D(0, 0);
            for (int i = 0; i <= degree_; i++)
            {
                Point2D control_point = new Point2D(pts_[i].x, pts_[i].y);
                sum.x += control_point.x * PascalValues[degree_][i] * one_min_t_pows[degree_ - i] * t_pows[i];
                sum.y += control_point.y * PascalValues[degree_][i] * one_min_t_pows[degree_ - i] * t_pows[i];
            }

            return sum;
        }

        private Point2D Bezier_DeCastlejau(float t)
        {
            List<Point2D> decast_vals = new List<Point2D>();
            float one_min_t = 1.0f - t;

            for (int i = 0; i < pts_.Count; ++i)
            {
                decast_vals.Add(pts_[i]);
            }

            for (int deg = degree_; deg >= 0; --deg)
            {
                for (int i = 0; i < deg; i++)
                {
                    decast_vals[i] = decast_vals[i] * one_min_t + decast_vals[i + 1] * t;
                }
            }

            return decast_vals[0];
        }



        private const float MAX_DIST = 6.0F;

        private void DrawBezier_Bern(System.Drawing.Graphics gfx, float alpha)
        {
            Point2D current_left;
            Point2D current_right = new Point2D(Bezier_Bernstein(0));

            for (float t = alpha; t < 1; t += alpha)
            {
                current_left = current_right;
                current_right = Bezier_Bernstein(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), Bezier_Bernstein(1).P());
        }

        private void DrawBezier_DeCast(System.Drawing.Graphics gfx, float alpha)
        {
            Point2D current_left;
            Point2D current_right = new Point2D(Bezier_DeCastlejau(0));

            for (float t = alpha; t < 1; t += alpha)
            {
                current_left = current_right;
                current_right = Bezier_DeCastlejau(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), Bezier_DeCastlejau(1).P());
        }

        private void DrawCircle(System.Drawing.Graphics gfx, System.Drawing.Pen pen, Point2D center, float radius=5.0f)
        {
            float x = (float)center.x - radius;
            float y = (float)center.y - radius;
            float width = 2 * radius;
            float height = 2 * radius;
            gfx.DrawEllipse(pen, x, y, width, height);
        }

        private const int maxRecv = 4;
        private void DrawMidpoint(System.Drawing.Graphics gfx, List<Point2D> cPs, int recvStep)
        {
            if (recvStep == maxRecv)
            {
                Point2D current_left;
                Point2D current_right = cPs[0];

                for (int i = 0; i < cPs.Count; ++i)
                {
                    current_left = current_right;
                    current_right = cPs[i];
                    gfx.DrawLine(splinePen, current_left.P(), current_right.P());
                }

                return;
            }

            List<Point2D> decast_vals = new List<Point2D>();
            const float one_half = 0.5f;

            List<Point2D> leftPoints = new List<Point2D>();
            List<Point2D> rightPoints = new List<Point2D>();

            for (int i = 0; i < cPs.Count; ++i)
            {
                decast_vals.Add(cPs[i]);
            }
            for (int deg = cPs.Count-1; deg >= 0; --deg)
            {
                leftPoints.Add(decast_vals[0]);
                for (int i = 0; i < deg; ++i)
                {
                    decast_vals[i] = decast_vals[i] * one_half + decast_vals[i + 1] * one_half;
                }
                rightPoints.Add(decast_vals[deg]);
            }

            DrawMidpoint(gfx, leftPoints, recvStep + 1);
            DrawMidpoint(gfx, rightPoints, recvStep + 1);
        }

        private void DrawDeCastlejau(System.Drawing.Graphics gfx, float alpha)
        {
            PointF p1 = new PointF(50, 50);
            PointF p2 = new PointF(950, 50);
            gfx.DrawLine(shellPen, p1, p2);

            p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, 3)).P();
            p2 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, -3)).P();

            gfx.DrawLine(yPen, p1, p2);

            p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, 0)).P();
            p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, 0)).P();
            gfx.DrawLine(xPen, p1, p2);
            gfx.DrawString("0", bFont, Brushes.Red, p1.X - 20, p1.Y - 10);

            for (int i = 1; i <= 3; i++)
            {
                p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, i)).P();
                p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, i)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                gfx.DrawString(i.ToString(), bFont, Brushes.Black, p1.X - 20, p1.Y - 10);

                p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, -i)).P();
                p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, -i)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                gfx.DrawString((-i).ToString(), bFont, Brushes.Black, p1.X - 20, p1.Y - 10);
            }
            foreach (var point in pts_)
            {
                float t = (float)WindowPointToNonNegativeGraphPoint(point).y;
                gfx.DrawString(t.ToString("F"), bFont, Brushes.Gray, (float)point.x, (float)point.y + 10);
            }

            Point2D current_left;
            Point2D current_right = new Point2D(DeCastlejau(0));

            for (float t = alpha; t < 1; t += alpha)
            {
                current_left = current_right;
                current_right = DeCastlejau(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), DeCastlejau(1).P());
        }

        private void DrawBernstein(System.Drawing.Graphics gfx, float alpha)
        {
            PointF p1 = new PointF(50, 50);
            PointF p2 = new PointF(950, 50);
            gfx.DrawLine(shellPen, p1, p2);

            p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, 3)).P();
            p2 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, -3)).P();

            gfx.DrawLine(yPen, p1, p2);

            p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, 0)).P();
            p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, 0)).P();
            gfx.DrawLine(xPen, p1, p2);
            gfx.DrawString("0", bFont, Brushes.Red, p1.X - 20, p1.Y - 10);

            for (int i = 1; i <= 3; i++)
            {
                p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, i)).P();
                p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, i)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                gfx.DrawString(i.ToString(), bFont, Brushes.Black, p1.X - 20, p1.Y - 10);

                p1 = NonNegativeGraphPointToWindowPoint(new Point2D(0.0f, -i)).P();
                p2 = NonNegativeGraphPointToWindowPoint(new Point2D(1.0f, -i)).P();
                gfx.DrawLine(hlinePen, p1, p2);
                gfx.DrawString((-i).ToString(), bFont, Brushes.Black, p1.X - 20, p1.Y - 10);
            }

            foreach (var point in pts_)
            {
                float t = (float)WindowPointToNonNegativeGraphPoint(point).y;
                gfx.DrawString(t.ToString("F"), bFont, Brushes.Gray, (float)point.x, (float)point.y + 10);
            }

            Point2D current_left;
            Point2D current_right = new Point2D(Bernstein(0));

            for (float t = alpha; t < 1; t += alpha)
            {
                current_left = current_right;
                current_right = Bernstein(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), Bernstein(1).P());
        }

        private void DrawPolyInterpolate(System.Drawing.Graphics gfx, float alpha)
        {
            gfx.DrawString("Interpolating Polynomial", bFont, Brushes.Black, 10, 40);

            if (pts_.Count < 2)
                return;
            ReComputeNewton();

            Point2D current_left;
            Point2D current_right = new Point2D(PolyInterpolate(0));

            for (float t = alpha; t < pts_.Count - 1; t += alpha)
            {
                current_left = current_right;
                current_right = PolyInterpolate(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), PolyInterpolate(pts_.Count - 1).P());
        }
        private void ReComputeNewton()
        {
            NewtonFormList.Clear();
            NewtonFormList.Add(pts_);
            for (int i = 1; i < pts_.Count; i++)
            {
                NewtonFormList.Add(new List<Point2D>());
                for (int j = 0; j < NewtonFormList[i-1].Count-1; j++)
                {
                    NewtonFormList[i].Add((NewtonFormList[i - 1][j + 1] - NewtonFormList[i - 1][j]) / (double)i);
                }
            }
        }
        private Point2D PolyInterpolate(float t)
        {
            Point2D result = new Point2D(0, 0);
            double t_val = 1.0;
            for (int i = 0; i < NewtonFormList.Count; i++)
            {
                result += NewtonFormList[i][0] * t_val;
                t_val *= (t - i);
            }

            return result;
        }
        private int PerformOperation(List<List<double>> a, int n)
        {
            int i, j, k = 0, c, flag = 0;

            // Performing elementary operations
            for (i = 0; i < n; i++)
            {
                if (a[i][i] == 0)
                {
                    c = 1;
                    while ((i + c) < n && a[i + c][i] == 0)
                        c++;
                    if ((i + c) == n)
                    {
                        flag = 1;
                        break;
                    }
                    for (j = i, k = 0; k <= n; k++)
                    {
                        double temp = a[j][k];
                        a[j][k] = a[j + c][k];
                        a[j + c][k] = temp;
                    }
                }

                for (j = 0; j < n; j++)
                {

                    // Excluding all i == j
                    if (i != j)
                    {

                        // Converting Matrix to reduced row
                        // echelon form(diagonal matrix)
                        double p = a[j][i] / a[i][i];

                        for (k = 0; k <= n; k++)
                            a[j][k] = a[j][k] - (a[i][k]) * p;
                    }
                }
            }
            return flag;
        }

        void GetResult(List<List<double>> a, List<double> coeff)
        {
            int n = a.Count;
            for (int i = 0; i < n; i++)
                coeff.Add(a[i][n] / a[i][i]);
        }
        private void ComputeSplineSolution()
        {
            List<List<double>> matrixX = new List<List<double>>();
            List<List<double>> matrixY = new List<List<double>>();

            List<(double, double)> tValues = new List<(double, double)>();
            for (int i = 0; i < pts_.Count; i++)
            {
                double t_2 = i * i;
                tValues.Add((t_2, t_2 * i));
            }

            int k_last = pts_.Count - 1;

            for (int i = 0; i < pts_.Count; i++)
            {
                //a0 a1 a2 a3
                (double t_2, double t_3) = tValues[i];
                List<double> rowX = new List<double> { 1.0, i, t_2, t_3 };
                List<double> rowY = new List<double> { 1.0, i, t_2, t_3 };
                for (int k = 1; k < k_last; k++)
                {
                    double piecewise_value = 0.0f;
                    if (i > k)
                    {
                        double t_min_c = i - k;
                        piecewise_value = t_min_c * t_min_c * t_min_c;
                    }

                    rowX.Add(piecewise_value);
                    rowY.Add(piecewise_value);
                }

                rowX.Add(pts_[i].x);
                rowY.Add(pts_[i].y);

                matrixX.Add(rowX);
                matrixY.Add(rowY);
            }

            List<double> derivative_x_zero = new List<double> { 0, 0, 2, 0 };
            List<double> derivative_y_zero = new List<double> { 0, 0, 2, 0 };
            for (int k = 1; k < k_last; k++)
            {
                derivative_x_zero.Add(0);
                derivative_y_zero.Add(0);
            }
            derivative_x_zero.Add(0);
            derivative_y_zero.Add(0);

            matrixX.Add(derivative_x_zero);
            matrixY.Add(derivative_y_zero);

            List<double> derivative_x_k = new List<double> { 0, 0, 2, 6 * k_last };
            List<double> derivative_y_k = new List<double> { 0, 0, 2, 6 * k_last };

            for (int k = 1; k < k_last; k++)
            {
                double innerVal = 6.0 * (k_last - k);
                derivative_x_k.Add(innerVal);
                derivative_y_k.Add(innerVal);
            }

            derivative_x_k.Add(0);
            derivative_y_k.Add(0);

            matrixX.Add(derivative_x_k);
            matrixY.Add(derivative_y_k);

            PerformOperation(matrixX, matrixX.Count);
            PerformOperation(matrixY, matrixY.Count);

            splineXCoeff.Clear();
            splineYCoeff.Clear();

            GetResult(matrixX, splineXCoeff);
            GetResult(matrixY, splineYCoeff);
        }
        private void DrawSplineInterpolate(System.Drawing.Graphics gfx, float alpha)
        {
            if (pts_.Count < 2)
                return;
            ComputeSplineSolution();

            Point2D current_left;
            Point2D current_right = new Point2D(SplineInterpolate(0));
            
            int k = pts_.Count - 1;

            for (float t = alpha; t < k; t += alpha)
            {
                current_left = current_right;
                current_right = SplineInterpolate(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }

            gfx.DrawLine(splinePen, current_right.P(), SplineInterpolate(k).P());
        }

        private Point2D SplineInterpolate(float t)
        {

            int k = pts_.Count - 1;
            double t_two = t * t;
            double t_three = t_two * t;
            double x = splineXCoeff[0] + splineXCoeff[1] * t + splineXCoeff[2] * t_two + splineXCoeff[3] * t_three;
            double y = splineYCoeff[0] + splineYCoeff[1] * t + splineYCoeff[2] * t_two + splineYCoeff[3] * t_three;

            for (int i = 1; i < k; i++)
            {
                double computed_value = 0;
                if (t >= i)
                    computed_value = t - i;
                x += splineXCoeff[3 + i] * computed_value * computed_value * computed_value;
                y += splineYCoeff[3 + i] * computed_value * computed_value * computed_value;

            }

            return new Point2D(x, y);
        }

        private void DrawDeBoor(System.Drawing.Graphics gfx, float alpha)
        {

            foreach (var point in pts_)
            {
                double x = WindowPointToGraphPoint(point).x;
                double y = WindowPointToGraphPoint(point).y;
                gfx.DrawString($"({x:F},{y:F})", Font, Brushes.Gray, (float)point.x, (float)point.y + 10);
            }

            Point2D current_left;
            Point2D current_right = new Point2D(DeBoorAlgthm(knot_[degree_]));
            double lastT = knot_[knot_.Count - degree_ - 1] - 0.01;
            for (float t = knot_[degree_]; t < lastT-0.1; t += alpha)
            {
                current_left = current_right;
                current_right = DeBoorAlgthm(t);
                gfx.DrawLine(splinePen, current_left.P(), current_right.P());
            }
            gfx.DrawLine(splinePen, current_right.P(), DeBoorAlgthm(lastT).P());
        }
        private Point2D DeBoorAlgthm(double t)
        {
            int s = pts_.Count - 1;
            int j = FindJ(t);
            List<Point2D> dynamicTable = new List<Point2D>(pts_);
            
            for (int k = 1; k < degree_; k++)
            {
                for (int i = s; i >= k; --i) 
                {
                    double coeffP_i = (t - knot_[i]) / (knot_[i + degree_ - (k - 1)] - knot_[i]);
                    double coeffP_i_1 = (knot_[i + degree_ - (k - 1)] - t) / (knot_[i + degree_ - (k - 1)] - knot_[i]);

                    dynamicTable[i] = dynamicTable[i] * coeffP_i + dynamicTable[i - 1] * coeffP_i_1;
                }
            }

            double coeffP_j = (t - knot_[j]) / (knot_[j + degree_ - (degree_ - 1)] - knot_[j]);
            double coeffP_j_1 = (knot_[j + degree_ - (degree_ - 1)] - t) / (knot_[j + degree_ - (degree_ - 1)] - knot_[j]);
            return dynamicTable[j] * coeffP_j + dynamicTable[j - 1] * coeffP_j_1;
        }

        private Point2D DeBoorAlgthmHelper(List<List<Point2D>> dynamicTable, List<List<bool>> isComputed, int k, int i, double t)
        {
            /*if (k == 0)
                return pts_[i];*/
            if (isComputed[k][i])
                return dynamicTable[k][i];

            double coeffP_i = (t - knot_[i]) / (knot_[i + degree_ - (k - 1)] - knot_[i]);
            double coeffP_i_1 = (knot_[i + degree_ - (k - 1)] - t) / (knot_[i + degree_ - (k - 1)] - knot_[i]);

            Point2D P_i_1 = DeBoorAlgthmHelper(dynamicTable, isComputed, k - 1, i - 1, t);
            Point2D P_i = DeBoorAlgthmHelper(dynamicTable, isComputed, k - 1, i, t);

            dynamicTable[k][i] = P_i * coeffP_i + P_i_1 * coeffP_i_1;
            return dynamicTable[k][i];

        }

        private int FindJ(double t)
        {
            if (t >= knot_[knot_.Count - 1])
                return knot_.Count - 1;

            int j = 0;
            
            for (int i = 0; i < knot_.Count; i++)
            {
                if (knot_[i] <= t && t < knot_[i + 1])
                {
                    j = i;
                    break;
                }
            }
            return j;
        }

        private void MAT290_Load(object sender, EventArgs e)
        {

        }

        private void SavePoints_Click(object sender, EventArgs e)
        {
            SaveFileDialog saveFileDialog = new SaveFileDialog();
            saveFileDialog.Filter = "Text|*.txt";
            saveFileDialog.Title = "Save an Points File";
            saveFileDialog.ShowDialog();

            if (saveFileDialog.FileName != "")
            {
                System.IO.FileStream fs =
                    (System.IO.FileStream)saveFileDialog.OpenFile();

                using (var writer=new StreamWriter(fs))
                {
                    writer.WriteLine($"{degree_}");

                    writer.WriteLine($"p");
                    foreach (var point in pts_)
                    {
                        writer.WriteLine($"{point.x} {point.y}");
                    }

                    writer.WriteLine($"k");
                    if (Menu_DeBoor.Checked)
                    {
                        foreach (var knot in knot_)
                        {
                            writer.WriteLine($"{knot}");
                        }
                    }
                }

                fs.Close();
            }

        }

        private void LoadPoints_Click(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog()
            {
                FileName = "Select a point file",
                Filter = "point files (*.txt)|*.txt",
                Title = "Open your saved points file"
            };

            if (openFileDialog.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    var filePath = openFileDialog.FileName;
                    using (FileStream fs = (FileStream)openFileDialog.OpenFile()) 
                    using (StreamReader sr = new StreamReader(fs))
                    {
                        //Process.Start("notepad.exe", filePath);
                        
                        var newpts = new List<Point2D>();
                        var newknot = new List<float>();
                        string line = null;
                        line = sr.ReadLine();
                        degree_ = Convert.ToInt32(line);
                        bool saveKnot = false;
                        while (true)
                        {
                            line=sr.ReadLine();
                            if(line==null)
                                break;
                            if (line == "p")
                            {
                                continue;
                            }
                            else if (line == "k")
                            {
                                saveKnot = true;
                                continue;
                            }

                            if (!saveKnot)
                            {
                                string[] words = line.Split(' ');
                                var x = Convert.ToDouble(words[0]);
                                var y = Convert.ToDouble(words[1]);
                                newpts.Add(new Point2D(x, y));
                            }
                            else
                            {
                                newknot.Add(Convert.ToSingle(line));
                            }
                            
                            //Console.WriteLine(line);
                        }

                        if (newpts.Count != 0)
                            pts_ = newpts;

                        if (newknot.Count != 0)
                            knot_ = newknot;


                        ResetKnotSeq();
                        UpdateKnotSeq();
                        NUD_degree.Value = degree_;
                        Refresh();
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show($"Error happend: check your file");
                }
            }
        }
    }
}