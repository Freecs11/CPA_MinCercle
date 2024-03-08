package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;

import supportGUI.Circle;
import supportGUI.Line;

public class DefaultTeam {
  public Line calculDiametre(ArrayList<Point> points) {
    // return tme1exercice7(points);
    return tme1exercice6(points);
}

// calculCercleMin: ArrayList<Point> --> Circle
//   renvoie un cercle couvrant tout point de la liste, de rayon minimum.
public Circle calculCercleMin(ArrayList<Point> points) {
    return tme1exercice5(points);
    //return tme1exercice4(points);
}
private Circle tme1exercice4(ArrayList<Point> inputPoints){
    ArrayList<Point> points = (ArrayList<Point>) inputPoints.clone();
    if (points.size()<1) return null;
    double cX,cY,cRadius,cRadiusSquared;
    for (Point p: points){
        for (Point q: points){
            cX = .5*(p.x+q.x);
            cY = .5*(p.y+q.y);
            cRadiusSquared = 0.25*((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
            boolean allHit = true;
            for (Point s: points)
                if ((s.x-cX)*(s.x-cX)+(s.y-cY)*(s.y-cY)>cRadiusSquared){
                    allHit = false;
                    break;
                }
            if (allHit) return new Circle(new Point((int)cX,(int)cY),(int)Math.sqrt(cRadiusSquared));
        }
    }
    double resX=0;
    double resY=0;
    double resRadiusSquared=Double.MAX_VALUE;
    for (int i=0;i<points.size();i++){
        for (int j=i+1;j<points.size();j++){
            for (int k=j+1;k<points.size();k++){
                Point p=points.get(i);
                Point q=points.get(j);
                Point r=points.get(k);
                //si les trois sont colineaires on passe
                if ((q.x-p.x)*(r.y-p.y)-(q.y-p.y)*(r.x-p.x)==0) continue;
                //si p et q sont sur la meme ligne, ou p et r sont sur la meme ligne, on les echange
                if ((p.y==q.y)||(p.y==r.y)) {
                    if (p.y==q.y){
                        p=points.get(k); //ici on est certain que p n'est sur la meme ligne de ni q ni r
                        r=points.get(i); //parce que les trois points sont non-colineaires
                    } else {
                        p=points.get(j); //ici on est certain que p n'est sur la meme ligne de ni q ni r
                        q=points.get(i); //parce que les trois points sont non-colineaires
                    }
                }
                //on cherche les coordonnees du cercle circonscrit du triangle pqr
                //soit m=(p+q)/2 et n=(p+r)/2
                double mX=.5*(p.x+q.x);
                double mY=.5*(p.y+q.y);
                double nX=.5*(p.x+r.x);
                double nY=.5*(p.y+r.y);
                //soit y=alpha1*x+beta1 l'equation de la droite passant par m et perpendiculaire a la droite (pq)
                //soit y=alpha2*x+beta2 l'equation de la droite passant par n et perpendiculaire a la droite (pr)
                double alpha1=(q.x-p.x)/(double)(p.y-q.y);
                double beta1=mY-alpha1*mX;
                double alpha2=(r.x-p.x)/(double)(p.y-r.y);
                double beta2=nY-alpha2*nX;
                //le centre c du cercle est alors le point d'intersection des deux droites ci-dessus
                cX=(beta2-beta1)/(double)(alpha1-alpha2);
                cY=alpha1*cX+beta1;
                cRadiusSquared=(p.x-cX)*(p.x-cX)+(p.y-cY)*(p.y-cY);
                if (cRadiusSquared>=resRadiusSquared) continue;
                boolean allHit = true;
                for (Point s: points)
                    if ((s.x-cX)*(s.x-cX)+(s.y-cY)*(s.y-cY)>cRadiusSquared){
                        allHit = false;
                        break;
                    }
                if (allHit) {System.out.println("Found r="+Math.sqrt(cRadiusSquared));resX=cX;resY=cY;resRadiusSquared=cRadiusSquared;}
            }
        }
    }
    return new Circle(new Point((int)resX,(int)resY),(int)Math.sqrt(resRadiusSquared));
}
private Circle tme1exercice5(ArrayList<Point> points){
    if (points.size()<1) return null;
    ArrayList<Point> rest = (ArrayList<Point>)points.clone();
    Point dummy=rest.get(0);
    Point p=dummy;
    for (Point s: rest) if (dummy.distance(s)>dummy.distance(p)) p=s;
    Point q=p;
    for (Point s: rest) if (p.distance(s)>p.distance(q)) q=s;
    double cX=.5*(p.x+q.x);
    double cY=.5*(p.y+q.y);
    double cRadius=.5*p.distance(q);
    rest.remove(p);
    rest.remove(q);
    while (!rest.isEmpty()){
        Point s=rest.remove(0);
        double distanceFromCToS=Math.sqrt((s.x-cX)*(s.x-cX)+(s.y-cY)*(s.y-cY));
        if (distanceFromCToS<=cRadius) continue;
        double cPrimeRadius=.5*(cRadius+distanceFromCToS);
        double alpha=cPrimeRadius/(double)(distanceFromCToS);
        double beta=(distanceFromCToS-cPrimeRadius)/(double)(distanceFromCToS);
        double cPrimeX=alpha*cX+beta*s.x;
        double cPrimeY=alpha*cY+beta*s.y;
        cRadius=cPrimeRadius;
        cX=cPrimeX;
        cY=cPrimeY;
    }
    return new Circle(new Point((int)cX,(int)cY),(int)cRadius);
}
private Line tme1exercice6(ArrayList<Point> points) {
    if (points.size()<2) return null;
    Point p=points.get(0);
    Point q=points.get(1);
    for (Point s: points) for (Point t: points) if (s.distance(t)>p.distance(q)) {p=s;q=t;}
    return new Line(p,q);
}

  public double produitVectorial(Point p1, Point p2, Point p3) {
    return (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x);
  }

  public double produitScalaire(Point p1, Point p2, Point p3) {
    return (p2.x-p1.x)*(p3.x-p1.x)+(p2.y-p1.y)*(p3.y-p1.y);
  }

  public ArrayList<Point> algo_Jarvis(ArrayList<Point> points){
    if (points.size()<3) {
      return null;
    }

    ArrayList<Point> enveloppe = new ArrayList<Point>();
    Point pointP = points.stream().min(Comparator.comparing(Point::getX)).get();
    Point firstPoint = pointP;
    Point virtuePoint = new Point(pointP.x, pointP.y-1);
    Point pointQ = null;
    for (Point point : points) {
      if (pointQ == null ){
        pointQ = point;
      }
      if (point.y < produitScalaire(pointP, virtuePoint, point)) {
        virtuePoint = point;
        pointQ = point;
      }
    }
    Point pointR = null;
    do{
      
      Point pointMin = null;
      double angleMin = Double.MAX_VALUE;
      for (Point p1 : points){
        if (p1 != pointP && p1 != pointQ) {
          double angle = Math.atan2(p1.y-pointP.y, p1.x-pointP.x) - Math.atan2(pointQ.y-pointP.y, pointQ.x-pointP.x);
          if (angle < 0) {
            angle += 2*Math.PI;
          }
          if (angle < angleMin) {
            angleMin = angle;
            pointMin = p1;
          }
        }
      }
      pointR = pointMin;
      pointP = pointQ;
      pointQ = pointR;
      
    }
    while (pointR != firstPoint);



    enveloppe.add(pointP);
    

    return enveloppe;
  }

  public ArrayList<Point> tri_pixel(ArrayList<Point> points){
    Point maxPoint = points.stream().max(Comparator.comparing(Point::getX)).get();

    Point[] pointsTriesYmin = new Point[maxPoint.x+1];
    for (int i = 0; i < points.size(); i++) {
      Point point = points.get(i);
      if (pointsTriesYmin[point.x] == null ){
        pointsTriesYmin[point.x] = point;
      } else {
        if (point.y < pointsTriesYmin[point.x].y) {
          pointsTriesYmin[point.x] = point;
        }
      }
    }

    Point[] pointsTriesYmax = new Point[maxPoint.x+1];
    for (int i = 0; i < points.size(); i++) {
      Point point = points.get(i);
      if (pointsTriesYmax[point.x] == null) {
        pointsTriesYmax[point.x] = point;
      } else {
        if (point.y > pointsTriesYmax[point.x].y) {
          pointsTriesYmax[point.x] = point;
        }
      }
    }

    ArrayList<Point> pointsTriesminY = new ArrayList<Point>();
    for (int i = 0; i < pointsTriesYmin.length; i++) {
      if (pointsTriesYmin[i] != null) {
        pointsTriesminY.add(pointsTriesYmin[i]);
      }
    }

    ArrayList<Point> pointsTriesmaxY = new ArrayList<Point>();
    for (int i = 0; i < pointsTriesYmax.length; i++) {
      if (pointsTriesYmax[i] != null) {
        pointsTriesmaxY.add(pointsTriesYmax[i]);
      }
    }

    Set<Point> pointsTries = new HashSet<>();
    pointsTries.addAll(pointsTriesminY);
    pointsTries.addAll(pointsTriesmaxY);
    ArrayList<Point> Pointstest  = new ArrayList<>(pointsTries);
    return Pointstest;
  }

  // enveloppeConvexe: ArrayList<Point> --> ArrayList<Point>
  //   renvoie l'enveloppe convexe de la liste.
  public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points){
    if (points.size()<3) {
      return null;
    }

    ArrayList<Point> enveloppe = new ArrayList<Point>();
    /*******************
     * PARTIE A ECRIRE *
     *******************/
    // Trie par pixel 
    ArrayList<Point> Pointstest = algo_Jarvis(points);
    return Pointstest;
  }

  public ArrayList<Point> algoNaif(ArrayList<Point> points){
    if (points.size()<3) {
      return null;
    }

    ArrayList<Point> enveloppe = new ArrayList<Point>();
    for (int i = 0; i < points.size(); i++) {
      for (int j = 0; j < points.size(); j++) {
        boolean flag = true;
        int countPositive = 0;
        int countNegative = 0;
        Point pointP = points.get(i);
        Point pointQ = points.get(j);
        for (int k = 0; k < points.size(); k++) {
          if (k != i && k != j) {

            Point pointR = points.get(k);
            if (produitVectorial(pointP, pointQ, pointR) > 0) {
              countPositive++;
            } else if (produitVectorial(pointP, pointQ, pointR) < 0) {
              countNegative++;
            } else {
              double distancePQ = pointP.distance(pointQ);
              double distancePR = pointP.distance(pointR);
              double distanceQR = pointQ.distance(pointR);
              if (distancePQ > distancePR && distancePQ > distanceQR) {
                continue;
              } else if (distancePR > distancePQ && distancePR > distanceQR) {
                pointQ = pointR;
              } else if (distanceQR > distancePQ && distanceQR > distancePR) {
                pointP = pointR;
              } else if (distancePQ == distancePR && distancePQ > distanceQR) {
                pointQ = pointR;
              } else if (distancePQ == distanceQR && distancePQ > distancePR) {
                pointQ = pointR;
              } else if (distancePR == distanceQR && distancePR > distancePQ) {
                pointP = pointR;
              } else if (distancePQ == distancePR && distancePQ == distanceQR) {
                pointQ = pointR;
              } else {
                continue;
              }
            }
            if (countNegative > 0 && countPositive > 0) {
              flag = false;
              break;
            }
          }
        }
        if (flag) {
          enveloppe.add(pointP);
          enveloppe.add(pointQ);
        }
      }
    } 

    
    
    return enveloppe;
  }



}
