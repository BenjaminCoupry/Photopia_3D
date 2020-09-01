using System;
using System.Collections.Generic;
using System.Drawing;

namespace Photopia_3D
{
    delegate Vector Mapping(double x, double y);
    delegate Vector AttenuationAngulaire(double teta);
    delegate EnergieLumineuse ConeEmission(double teta);

    delegate Vector FiltrationVolumetrique(Vector position);
    delegate Vector RepartitionAngulaireParametree(double teta, double FiltrationLocale);

    class Mapper
    {
        //public static Mapping fromBitmap(Bitmap bp)
        public static Mapping Uniforme(Vector valeur)
        {
            return (x, y) => valeur;
        }
        public static Mapping BPAngulaire(Bitmap valeurs)
        {
            return (x, y) => {
                double kteta = x / Math.PI;
                double kphy = y / Math.PI;
                Color c = valeurs.GetPixel((int)(kphy * (double)(valeurs.Width-1)), (int)(kteta * (double)(valeurs.Height-1)));
                return new Vector(c.R / 255.0, c.G / 255.0, c.B / 255.0);
            };
        }
        public static AttenuationAngulaire AttenuationGaussienne(double lambda1, double lambda2, double lambda3)
        {
            return (x) => new Vector(Math.Exp(-lambda1 * x), Math.Exp(-lambda2 * x), Math.Exp(-lambda3 * x));
        }
        public static AttenuationAngulaire Cos()
        {
            return (x) => new Vector(1, 1, 1) * Math.Max(0, Math.Cos(x));
        }
        public static AttenuationAngulaire Sin()
        {
            return (x) => new Vector(1,1,1)* Math.Max(0, Math.Sin(x));
        }
        public static ConeEmission ConeCosinus(Color Base)
        {
            EnergieLumineuse E0 = new EnergieLumineuse(Base);
            return (x) => E0*Math.Max(0,Math.Cos(x));
        }
        public static ConeEmission ConeConstant(Color Base)
        {
            EnergieLumineuse E0 = new EnergieLumineuse(Base);
            return (x) => E0 ;
        }
        public static FiltrationVolumetrique FiltrationUniforme(Vector Filtre)
        {
            return (x) => Filtre;
        }
        public static RepartitionAngulaireParametree Cacahuete(double Coeff, double A0)
        {
            
            return (x,y) => Math.Pow(Math.Abs(Math.Cos(x)),Coeff/Math.Abs(y))*A0*new Vector(1,1,1);
        }
        public static RepartitionAngulaireParametree Nulle()
        {

            return (x, y) => new Vector(0,0,0);
        }

    }

    class Materiau
    {
        public Mapping Coloration;
        public Mapping Normales;
        public AttenuationAngulaire AttenuationReflexion;
        public AttenuationAngulaire ColorationAngulaire;
        public Vector FiltrationVolumetrique;
        public double n_refraction;
        public bool Opaque;

        public Materiau(Mapping coloration, Mapping normales, AttenuationAngulaire attenuationReflexion, AttenuationAngulaire colorationAngulaire, Vector filtrationVolumetrique, double n_refraction, bool opaque)
        {
            Coloration = coloration;
            Normales = normales;
            AttenuationReflexion = attenuationReflexion;
            ColorationAngulaire = colorationAngulaire;
            FiltrationVolumetrique = filtrationVolumetrique;
            this.n_refraction = n_refraction;
            Opaque = opaque;
        }

        public Vector NormaleReelle(Tuple<double, double> coord2D, Base baseLocale, Vector NormaleGeometrique)
        {
            return (NormaleGeometrique + (baseLocale.Delocaliser(Normales(coord2D.Item1, coord2D.Item2)) - baseLocale.position)).Normer();
        }
        public Vector CouleurReelle(Tuple<double, double> coord2D, InfoContact Contact, Vector Sortant)
        {
            double AngleSpeculaire = Contact.Rebond().direction.Angle(Sortant);
            double AngleIncidence = (-Contact.VecteurIncident).Angle(Contact.VecteurNormal);
            Vector CouleurLocale = Coloration(coord2D.Item1, coord2D.Item2);
            Vector FacteurRasant = ColorationAngulaire(AngleIncidence);
            if(FacteurRasant.isNan())
            {
                ;
            }
            if (CouleurLocale.isNan())
            {
                ;
            }
            Vector retour = AttenuationReflexion(AngleSpeculaire) % (new Vector(1, 1, 1) % FacteurRasant + CouleurLocale % (1.0 - FacteurRasant));
            if (retour.isNan())
            {
                ;
            }
            return retour;
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            Materiau Test = new Materiau(Mapper.Uniforme(new Vector(1, 0.1, 0.1)), Mapper.Uniforme(new Vector(0, 0, 0)), Mapper.AttenuationGaussienne(1, 1, 0.5), Mapper.AttenuationGaussienne(2, 2,2), new Vector(1, 1, 1), 1, true);
            Materiau Test2 = new Materiau(Mapper.Uniforme(new Vector(1, 0.1, 0.1)), Mapper.Uniforme(new Vector(0, 0, 0)), Mapper.AttenuationGaussienne(1, 1, 0.5), Mapper.AttenuationGaussienne(2, 2, 2), new Vector(0.5, 0.5, 0.5), 1, false);
            SourceLumineuse Source = new SourceLumineuse(new Base(Vector.X(), Vector.Y(), Vector.Z(), new Vector(0, 4, 3)), Mapper.ConeConstant(Color.White));
            SourceLumineuse Source2 = new SourceLumineuse(new Base(Vector.X(), Vector.Y(), Vector.Z(), new Vector(0, 4, -3)), Mapper.ConeConstant(Color.Red));
            Sphere S = new Sphere(1.8, Test, new Base(Vector.X(), Vector.Y(), Vector.Z(), new Vector(0, 0, 0)));
            Sphere S1 = new Sphere(0.5, Test2, new Base(Vector.X(), Vector.Y(), Vector.Z(), new Vector(3, 0, 0)));
            Univers U = new Univers(new List<ObjetPhysique>() { S,S1 }, new List<SourceLumineuse>() {  Source2,Source }, 3, 5, 0.01, 1, Mapper.FiltrationUniforme(new Vector(0.5, 0.7, 0.2)), Mapper.Cacahuete(0.5,0.2),0.2);
            Camera C = new Camera(new Base(Vector.X(), Vector.Y(), Vector.Z(), new Vector(0, 3, 0)), 2, 100, 100, 0.05);
            C.Rendre(U).Save("D:/lab/Rendu3D/1.bmp");
        }
    }
    class Vector
    {
        public double x;
        public double y;
        public double z;

        public Vector(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public double Norme()
        {
            return Math.Sqrt(x * x + y * y + z * z);
        }
        public static Vector operator +(Vector a, Vector b)
        {
            return new Vector(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        public static Vector operator +(double a, Vector b)
        {
            return new Vector(a + b.x, a + b.y, a + b.z);
        }
        public static Vector operator +(Vector a, double b)
        {
            return new Vector(a.x + b, a.y + b, a.z + b);
        }
        public static Vector operator -(double a, Vector b)
        {
            return new Vector(a - b.x, a - b.y, a - b.z);
        }
        public static Vector operator -(Vector a, double b)
        {
            return new Vector(a.x - b, a.y - b, a.z - b);
        }
        public static Vector operator -(Vector a, Vector b)
        {
            return new Vector(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        public static double operator *(Vector a, Vector b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
        public static Vector operator %(Vector a, Vector b)
        {
            return new Vector(a.x * b.x , a.y * b.y , a.z * b.z);
        }
        public static Vector operator -(Vector a)
        {
            return new Vector(-a.x, -a.y, -a.z);
        }
        public static Vector operator *(Vector a, double b)
        {
            return new Vector(a.x * b, a.y * b, a.z * b);
        }
        public static Vector operator *(double a, Vector b)
        {
            return b * a;
        }
        public static Vector operator /(Vector a, double b)
        {
            return a * (1.0 / b);
        }
        public static Vector operator ^(Vector a, Vector b)
        {
            return new Vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }
        public static Vector operator / (Vector a, Vector b)
        {
            return (a*(b.Normer()))*(b.Normer());
        }
        public static Vector Zero()
        {
            return new Vector(0, 0, 0);
        }
        public static Vector X()
        {
            return new Vector(1, 0, 0);
        }
        public static Vector Y()
        {
            return new Vector(0, 1, 0);
        }
        public static Vector Z()
        {
            return new Vector(0, 0, 1);
        }
        public Vector Normer()
        {
            if (this.Norme() != 0)
            {
                return this / this.Norme();
            }
            else
            {
                return Vector.Zero();
            }
        }
        public static double Distance(Vector a, Vector b)
        {
            return (a - b).Norme();
        }
        public double Angle(Vector b)
        {
            double n = this.Norme() * b.Norme();
            double retour;
            if (n != 0)
            {
                double Cos = Math.Min(1.0,Math.Max(-1.0,(this * b) / (n)));
                retour =  Math.Acos(Cos);
                if (double.IsNaN(retour))
                {
                    ;
                }
            }
            else
            {
                retour =  0;
            }
            
            return retour;
        }
        public bool isNan()
        {
            return double.IsNaN(x) || double.IsNaN(y) || double.IsNaN(z);
        }
    }
    class Localisable
    {
        public Base Orientation;
        public Vector Position()
        {
            return Orientation.position;
        }
        public Localisable( Base orientation)
        {
            Orientation = orientation;
        }
    }
    abstract class ObjetPhysique : Localisable
    {
        Materiau materiau;
        public ObjetPhysique(Materiau Mat, Base orientation) : base( orientation)
        {
            materiau = Mat;
        }
        public abstract double Distance(Vector Position);
        public abstract Tuple<double, double> getCoordonees2D(Vector position);
        public abstract bool Contact(Vector Position, double d0);
        public abstract Base GetBase(Tuple<double, double> coord2D);
        public abstract Vector Normale(Tuple<double, double> coord2D);
        public Vector GetCoordLocale(Vector position)
        {
            Base locale = Orientation;
            Vector coordLocales = locale.GetCoordoneesDansBase(position);
            return coordLocales;
        }
        public Vector NormaleReelle(Vector position)
        {
            Tuple<double, double> coord2D = getCoordonees2D(position);
            return getMateriau().NormaleReelle(coord2D,GetBase(coord2D),Normale(coord2D));
        }
        public Vector CouleurReelle( InfoContact Contact, Vector Sortant)
        {
            Vector retour = getMateriau().CouleurReelle(getCoordonees2D(Contact.Position), Contact, Sortant);
            if(retour.isNan())
            {
                ;
            }
            return retour;
        }
        public Materiau getMateriau()
        {
            return materiau;
        }
        public abstract Vector Repousser(Vector position, double dist);

    }
    class Base
    {
        public Vector u1;
        public Vector u2;
        public Vector u3;
        public Vector position;

        public Base(Vector u1, Vector u2, Vector u3, Vector position)
        {
            this.u1 = u1;
            this.u2 = u2;
            this.u3 = u3;
            this.position = position;
        }
        public Base()
        {
            this.u1 = Vector.X();
            this.u2 = Vector.Y();
            this.u3 = Vector.Z();
            this.position = Vector.Zero();
        }
        public Base(RepereSpherique sphere)
        {
            this.u1 = sphere.GetUphy();
            this.u2 = sphere.GetUteta();
            this.u3 = sphere.GetUr();
            this.position = sphere.GetCartesien();
        }

        public Vector Delocaliser(Vector coordoneesBase)
        {
            return position + u1 * coordoneesBase.x + u2 * coordoneesBase.y + u3 * coordoneesBase.z;
        }
        public Base Delocaliser(Base baseDel)
        {
            return new Base(Delocaliser(baseDel.u1) - position, Delocaliser(baseDel.u2) - position, Delocaliser(baseDel.u3) - position, Delocaliser(baseDel.position));
        }
        public Vector GetCoordoneesDansBase(Vector Cartesien)
        {
            Vector local = Cartesien - position;
            return new Vector(local * u1, local * u2, local * u3);
        }
    }
    class RepereSpherique
    {
        public double r;
        public double teta;
        public double phy;

        public RepereSpherique(double r, double teta, double phy)
        {
            this.r = r;
            this.teta = teta;
            this.phy = phy;
        }
        public RepereSpherique(Vector cartesien)
        {
            r = cartesien.Norme();
            teta = Math.Acos(cartesien.z / r);
            phy = Math.Atan2(cartesien.y, cartesien.x);
        }
        public Vector GetCartesien()
        {
            return new Vector(r * Math.Sin(teta) * Math.Cos(phy), r * Math.Sin(teta) * Math.Sin(phy), r * Math.Cos(teta));
        }
        public Vector GetUr()
        {
            return new Vector(Math.Sin(teta) * Math.Cos(phy), Math.Sin(teta) * Math.Sin(phy), Math.Cos(teta));
        }
        public Vector GetUteta()
        {
            return new Vector(Math.Cos(teta) * Math.Cos(phy), Math.Cos(teta) * Math.Sin(phy), -Math.Sin(teta));
        }
        public Vector GetUphy()
        {
            return new Vector(-Math.Sin(phy), Math.Cos(phy), 0);
        }
        public Base GetBase()
        {
            return new Base(GetUr(), GetUteta(), GetUphy(), GetCartesien());
        }
    }

    class Tir
    {
        public Vector direction;
        public Vector depart;
        public Tir(Vector direction, Vector depart)
        {
            this.direction = direction;
            this.depart = depart;
        }
    }
    class Rayon
    {
        public Tir tir;
        public double Longueur;
        public double Pas;
        public double distMax;
        public Rayon(Tir tir, double longueur, double Pas, double distMax)
        {
            this.tir = tir;
            Longueur = longueur;
            this.Pas = Pas;
            this.distMax = distMax;
        }
        public InfoContact Tirer(Univers univ, List<ObjetPhysique> Exclusion)
        {
            ObjetPhysique contact;
            do
            {
                Avancer(true,univ);
                contact = Contact(univ, Exclusion,univ.EpsilonContact);
            } while (contact is null && Longueur < distMax);
            if(contact is null)
            {
                return new InfoContact(Position(), tir.direction, new Vector(0,0,0), null);
            }
            
            return new InfoContact(Position(),tir.direction.Normer(),contact.NormaleReelle(Position()),contact);
        }
        public ObjetPhysique Contact(Univers univ, List<ObjetPhysique> Exclusion, double d0)
        {
            foreach(ObjetPhysique obj in univ.objets)
            {
                if(!Exclusion.Contains(obj))
                {
                    if(obj.Contact(this.Position(),d0))
                    {
                        return obj;
                    }
                }
            }
            return null;
        }
        Vector Opacite (Univers univ, double LongueurParcours)
        {
            //Calcule la somme d'opacite sur un parcours
            Vector retour = new Vector(0, 0, 0);
            ObjetPhysique contact;
            do
            {
                Avancer(false,univ);
                contact = Contact(univ, new List<ObjetPhysique>(),univ.EpsilonContact);
                if(contact is null)
                {
                    retour += univ.atmosphere(Position()) * Pas;
                }
                else
                {
                    if (contact.getMateriau().Opaque)
                    {
                        return new Vector(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
                    }
                    else
                    {
                        retour += contact.getMateriau().FiltrationVolumetrique * Pas;
                    }
                }
            } while (Longueur<LongueurParcours && Longueur < distMax);
            return retour;
        }
        public static Vector CalculerOpacite (Univers univ, Vector A, Vector B)
        {
            Rayon r = new Rayon(new Tir((B - A).Normer(), A),univ.PasVolumetrique, univ.PasVolumetrique, univ.distRayonMax);
            return r.Opacite(univ, (B - A).Norme()-2*univ.PasVolumetrique);
        }
        public static EnergieLumineuse CalculerAtmosphere(Univers univ, Tir tir)
        {
            Rayon r = new Rayon(tir, 0, univ.EpsilonContact, univ.distRayonMax);
            EnergieLumineuse E = new EnergieLumineuse();
            Vector A = tir.depart;
            while(r.Longueur<r.distMax)
            {
                Vector X = r.Position();
                foreach(SourceLumineuse sour in univ.sources)
                {
                    Vector B = sour.Position();
                    EnergieLumineuse recueBrute = sour.GetEnergieLumineuse(X);

                    Vector Opacite = Rayon.CalculerOpacite(univ, B, X) + Rayon.CalculerOpacite(univ, X, A);
                    EnergieLumineuse Opacifiee = recueBrute.Opacifier(Opacite, univ.facteurExtinction);
                    double angle = (X - B).Angle(A - X);
                    double densite = univ.atmosphere(X).Norme();
                    Vector FacteurReorientation = univ.reorientationAtmosphere(angle, densite);
                    E += (Opacifiee*FacteurReorientation)*r.Pas;

                }
                r.Avancer(false,univ);
            }
            return E;
        }
        public Vector Position()
        {
            return tir.depart + tir.direction * Longueur;
        }
        public void Avancer(bool vitesseMax,Univers univ)
        {
            //VitesseMax : avancer de la distance qui separe de l'objet le plus proche
            if (!vitesseMax)
            {
                Longueur += Pas;
            }
            else
            {
                double distmin = double.PositiveInfinity;
                foreach(ObjetPhysique obj in univ.objets)
                {
                    double dist = obj.Distance(Position());
                    if(dist<distmin)
                    {
                        distmin = dist;
                    }
                }
                Longueur += distmin;
            }
        }
    }
    class Camera : Localisable
    {

        double Recul;
        int X_ecran;
        int Y_ecran;
        double Scale;

        public Camera( Base Orientation, double recul, int x_ecran, int y_ecran, double scale) : base(Orientation)
        {
            Recul = recul;
            X_ecran = x_ecran;
            Y_ecran = y_ecran;
            Scale = scale;
        }

        public Tir GetTirPixel(int x, int y)
        {
            Vector PositionSurEcran = new Vector(x * Scale - X_ecran * Scale / 2.0,0, y * Scale - Y_ecran * Scale / 2.0);
            Vector PositionReelle = Orientation.Delocaliser(PositionSurEcran);
            Vector PositionPointFuite = Orientation.Delocaliser(Recul * Vector.Y());
            return new Tir((PositionReelle - PositionPointFuite).Normer(), PositionReelle);
        }

        public Bitmap Rendre(Univers univ)
        {
            Bitmap ret = new Bitmap(X_ecran, Y_ecran);
            for(int x=0;x<X_ecran;x++)
            {
                Console.WriteLine(x);
                for (int y = 0; y < Y_ecran; y++)
                {
                    Tir tir = GetTirPixel(x, y);
                    Rayon ray = new Rayon(tir, 0, univ.EpsilonContact,univ.distRayonMax);
                    EnergieLumineuse col = univ.GetLumiere(ray,univ.bounceLimit,new List<ObjetPhysique>());
                    ret.SetPixel(x, y, col.VersRGB());
                }
            }
            return ret;
        }
    }
    class SourceLumineuse : Localisable
    {
        public ConeEmission coneEmission;

        public SourceLumineuse(Base orientation, ConeEmission coneEmission_) : base(orientation)
        {
            coneEmission = coneEmission_;
        }
        public EnergieLumineuse GetEnergieLumineuse(Vector Position)
        {
            return coneEmission(Orientation.u1.Angle(Position - Orientation.position));
        }
    }
    class EnergieLumineuse
    {
        static double coeffCaptage = 1;
        public double r;
        public double g;
        public double b;
        public EnergieLumineuse()
        {
            r = 0;
            g = 0;
            b = 0;
        }
        public EnergieLumineuse(double r, double g, double b)
        {
            this.r = r;
            this.g = g;
            this.b = b;
        }
        public EnergieLumineuse(Color color)
        {
            r = convInv(color.R / 256.0);
            g = convInv(color.G / 256.0);
            b = convInv(color.B / 256.0);
        }
        public Color VersRGB()
        {
            return Color.FromArgb(Convert.ToInt32(255.0 * conv(r)), Convert.ToInt32(255.0 * conv(g)), Convert.ToInt32(255.0 * conv(b)));
        }
        private double conv(double x)
        {
            return 1.0 - Math.Exp(-coeffCaptage * x);
        }
        private double convInv(double x)
        {
            return -Math.Log(1.0-x)/coeffCaptage;
        }
        public static  EnergieLumineuse operator + (EnergieLumineuse a, EnergieLumineuse b)
        {
            return new EnergieLumineuse(a.r + b.r, a.g + b.g, a.b + b.b);
        }
        public static EnergieLumineuse operator *(EnergieLumineuse a, Vector b)
        {
            return new EnergieLumineuse(a.r * b.x, a.g * b.y, a.b * b.z);
        }
        public static EnergieLumineuse operator *(EnergieLumineuse a, double b)
        {
            return new EnergieLumineuse(a.r * b, a.g * b, a.b * b);
        }
        public EnergieLumineuse Opacifier(Vector Opacite, double facteurExtinction)
        {
            return new EnergieLumineuse(r * Math.Exp(-facteurExtinction * Opacite.x), g * Math.Exp(-facteurExtinction * Opacite.y), b * Math.Exp(-facteurExtinction * Opacite.z));
        }
        public bool isNan()
        {
            return double.IsNaN(r) || double.IsNaN(g) || double.IsNaN(b);
        }

    }
    class Univers
    {
        public List<ObjetPhysique> objets;
        public List<SourceLumineuse> sources;
        public int bounceLimit;
        public double distRayonMax;
        public double EpsilonContact;
        public double PasVolumetrique;
        public double facteurExtinction;
        public FiltrationVolumetrique atmosphere;
        public RepartitionAngulaireParametree reorientationAtmosphere;

        public Univers(List<ObjetPhysique> objets, List<SourceLumineuse> sources, int bounceLimit, double distRayonMax, double epsiloncontact, double facteurExtinction, FiltrationVolumetrique atmosphere, RepartitionAngulaireParametree reorientationAtmosphere, double pasVolumetrique)
        {
            this.objets = objets;
            this.sources = sources;
            this.bounceLimit = bounceLimit;
            this.distRayonMax = distRayonMax;
            EpsilonContact = epsiloncontact;
            this.facteurExtinction = facteurExtinction;
            this.atmosphere = atmosphere;
            this.reorientationAtmosphere = reorientationAtmosphere;
            PasVolumetrique = pasVolumetrique;
        }

        public EnergieLumineuse GetLumiere(Rayon R, int bounceLeft, List<ObjetPhysique> ExclusionContact)
        {
            EnergieLumineuse E = Rayon.CalculerAtmosphere(this, R.tir);
            if(E.isNan())
            {
                ;
            }
            EnergieLumineuse Ei = new EnergieLumineuse();
            InfoContact contact = R.Tirer(this, ExclusionContact);
            if(!(contact.Objet is null))
            {
                Vector Opacite = Rayon.CalculerOpacite(this, R.tir.depart, R.Position());
                Materiau mat = contact.Objet.getMateriau();
                if (bounceLeft > 0)
                {
                    Tir rebond = contact.Rebond();
                    Tir refrac = contact.Refraction();
                    Rayon Rebond = new Rayon(rebond, 0, R.Pas, R.distMax);
                    Rayon Refrac = new Rayon(refrac, 0, R.Pas, R.distMax);
                    EnergieLumineuse LumRebond = GetLumiere(Rebond, bounceLeft - 1, new List<ObjetPhysique>() { contact.Objet });
                    if (!mat.Opaque)
                    {
                        EnergieLumineuse LumRefrac = GetLumiere(Refrac, bounceLeft - 1, new List<ObjetPhysique>() { contact.Objet });
                        Ei += LumRefrac;
                    }
                    Vector Couleur = contact.Objet.CouleurReelle(contact, rebond.direction);
                    Ei += LumRebond * Couleur;

                }

                foreach (SourceLumineuse lu in sources)
                {
                    Vector OpaciteLum = Rayon.CalculerOpacite(this, lu.Position(), contact.Position);
                    Vector Sortant =  (lu.Position()- contact.Position).Normer();
                    EnergieLumineuse ArriveeBrute = lu.GetEnergieLumineuse(contact.Position);
                    EnergieLumineuse Attenuee = ArriveeBrute.Opacifier(OpaciteLum, facteurExtinction);
                    Vector Couleur = contact.Objet.CouleurReelle(contact,Sortant);
                    EnergieLumineuse PointContact = Attenuee * Couleur;
                    Ei += PointContact;
                }
                Ei.Opacifier(Opacite, facteurExtinction);
            }
            return E + Ei;
        }
    }

    class InfoContact
    {
        public Vector Position;
        public Vector VecteurIncident;
        public Vector VecteurNormal;
        public ObjetPhysique Objet;
        
        public bool Contact()
        {
            return !(Objet is null);
        }
        public Tir Rebond()
        {
            double u0 = VecteurNormal * (-VecteurIncident);
            Tir retour = new Tir(VecteurIncident+2*u0*VecteurNormal, Objet.Repousser(Position,0.01));
            return retour;
        }
        public Tir Refraction()
        {
            double n = Objet.getMateriau().n_refraction;
            double u0 = VecteurNormal * (-VecteurIncident);
            double u1 = 1.0 - Math.Pow(1.0 / n, 2) * (1.0 - u0 * u0);
            if(u1<0)
            {
                return Rebond();
            }
            else
            {
                Vector Retour;
                double u2 = Math.Sqrt(u1);
                if(u0>0)
                {
                    Retour = (1.0 / n) * VecteurIncident + ((1.0 / n) * u0 - u2) * VecteurNormal;
                }
                else
                {
                    Retour = (1.0 / n) * VecteurIncident + ((1.0 / n) * u0 + u2) * VecteurNormal;
                }
                return new Tir(Retour, Position);
            }
        }
        public InfoContact(Vector position, Vector vecteurIncident, Vector vecteurNormal, ObjetPhysique objet)
        {
            Position = position;
            VecteurIncident = vecteurIncident;
            VecteurNormal = vecteurNormal;
            Objet = objet;
        }
    }

    class Sphere : ObjetPhysique
    {
        double rayon;
        public Sphere(double Rayon, Materiau Mat, Base orientation) : base(Mat, orientation)
        {
            rayon = Rayon;
        }

        public override bool Contact(Vector Position, double d0)
        {
            return (Position-Orientation.position).Norme()<rayon+d0;
        }

        public override double Distance(Vector Position)
        {
            return Math.Max(0, (Orientation.position - Position).Norme() - rayon);
        }

        public override Base GetBase(Tuple<double, double> coord2D)
        {
            //(0,0,0) est le centre de la sphere
            RepereSpherique rep = new RepereSpherique(rayon,coord2D.Item1,coord2D.Item2);
            Base B = rep.GetBase();
            return B;
        }

        public override Tuple<double, double> getCoordonees2D(Vector position)
        {
            RepereSpherique rep = new RepereSpherique(GetCoordLocale(position));
            return new Tuple<double, double>(rep.teta, rep.phy);
        }

        public override Vector Normale(Tuple<double, double> coord2D)
        {
            RepereSpherique rep = new RepereSpherique(rayon, coord2D.Item1, coord2D.Item2);
            return Orientation.Delocaliser(rep.GetUr())-Orientation.position;
        }

        public override Vector Repousser(Vector position, double dist)
        {
            RepereSpherique rep = new RepereSpherique(GetCoordLocale(position));
            rep.r +=dist;
            return Orientation.Delocaliser(rep.GetCartesien());
        }
    }



}
