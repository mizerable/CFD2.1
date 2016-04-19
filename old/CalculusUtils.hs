{-# LANGUAGE ScopedTypeVariables, RankNTypes , KindSignatures , FlexibleContexts #-}
module CalculusUtils where

import SolutionDomain
import Data.List
         
direcDimenType:: Direction -> DimensionType
direcDimenType direc = case direc of
    Time -> Temporal
    _ -> Spatial
  
volumeOrInterval dimetype position vs = case dimetype of
    Temporal -> timeStep
    Spatial -> foldl' (\p d  -> p * sideLength d position vs) 1 $! enumFrom X 
    
volumeOrInterval_adj dimetype x env = case dimetype of
    Temporal -> timeStep
    Spatial -> foldl' (\p d  -> p * sideLength_adj d x env) 1 $! enumFrom X

c2D:: Property -> Direction
c2D c = case c of
    U -> X
    V -> Y
    W -> Z
    _ -> error "other properties don't have an associated direction"
    
d2C:: Direction -> Property
d2C d = case d of 
    X -> U
    Y -> V
    Z -> W
    _ -> error "other directions don't have an associated property"

--solveUnknown::(Fractional a)=> Equation (SchemeType ->Term a)->Position->a
solveUnknown (Equation l r _) position vs= 
    let sumUnknown p n =  p + case n of
            Unknown u-> u
            SubExpression s -> sumExpression sumUnknown $! getTerms s
            _ -> 0
        sumConstants p n =  p + case n of
            Constant c-> c
            --Derivative {}-> sumExpression sumConstants [approximateDerivative n position vs]  
            Derivative_adj {}-> sumExpression sumConstants [approximateDerivative_adj n position vs]
            SubExpression s -> sumExpression sumConstants $! getTerms s
            _ -> 0
        sumExpression s = foldl' s 0
        lhsUnknown = sumExpression sumUnknown l
        rhsUnknown = sumExpression sumUnknown r
        lhsConstants = sumExpression sumConstants l
        rhsConstants = sumExpression sumConstants r
    in (rhsConstants - lhsConstants)/(lhsUnknown-rhsUnknown)
    --in (1 - 0)/(1-0)

distributeMultiply::(Num a)=> [Term a]->a->[Term a]
distributeMultiply terms m = concatMap (multTerm m) terms
        
--multTerm:: (Num a)=> a -> Term a -> [Term a]
multTerm m term  = case term of
    SubExpression s -> distributeMultiply ( getTerms s  ) m
    {-Derivative direc func side multF->  
        let modf st x s = fmap (*m) (func st x s)   
        in [Derivative direc modf side multF]-}
    Derivative_adj direc func side multF->  
        let modf st x s = fmap (*m) (func st x s)   
        in [Derivative_adj direc modf side multF]
    _-> [fmap (*m) term]
 
integSurface_adj f position direction unknownProp vs= 
        let (s1, s2) = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f 
                        else (\x -> let [res] = multTerm (-1) x
                                in res  ).f 
                    term = modf s
                    sideAreaVal = sideArea_adj s position vs
                    isUnknown = (direcDimenType direction,isUpper) == (Temporal,True) 
                in case term of
                    Derivative_adj d subf _ m-> head $ multTerm sideAreaVal $ Derivative_adj d subf s m
                    SubExpression (Expression expr) -> SubExpression $ Expression $ distributeMultiply expr sideAreaVal
                    _-> if isUnknown 
                        then 
                            let constantVal = fmap
                                    (* (sideAreaVal 
                                        / prop_adj Directional unknownProp position s vs)) 
                                    term
                            in Unknown $ if isNaN (val constantVal) 
                                then 1 else val constantVal      
                        else fmap (* sideAreaVal) term
        in [value s1 True , value s2 False]       

integSingleTerm_adj term dimetype cellposition unknownProp vs=  
        let nonDerivAnswer = case term of 
                SubExpression _ -> error "can't integrate a subexpression as a single term"
                _ -> multTerm (volumeOrInterval_adj dimetype cellposition vs) term   
        in case term of     
            Derivative_adj direction func _ _->
                if direcDimenType direction == dimetype  
                then integSurface_adj ( func Directional cellposition ) cellposition direction unknownProp vs
                else nonDerivAnswer
            _ -> nonDerivAnswer    
         
integ_adj vs dimetype terms cellposition unknownProp= 
    case terms of
        [] -> []
        (x:xs) -> integSingleTerm_adj x dimetype cellposition unknownProp vs 
            ++ integ_adj vs dimetype xs cellposition unknownProp

d_::[Property]-> Direction -> AdjGraph  -> (Term Double)
d_ properties d= df_ properties 1 d 

-- df_ :: (Num a,Fractional a) => [Property]-> a ->Direction -> Reader (ValSet a) (Term a)
df_ :: [Property]-> Double ->Direction -> AdjGraph -> (Term Double)      
df_ properties factor d= dfm_ properties factor (\_ _ _ -> 1) d         

dfm_ properties factor multF direction d= 
    Derivative_adj direction 
        (\st x s -> Constant $  multProps_adj properties st x s d   * factor )
        Center
        multF      
      
--dd_ :: (Num a,Fractional a)=> Reader (ValSet a) (Term a) -> Direction -> Reader (ValSet a) (Term a)
dd_ inner  d= ddf_ inner 1 d 

--ddf_ ::(Num a,Fractional a)=> Reader (ValSet a) (Term a) -> a-> Direction -> Reader (ValSet a) (Term a)
ddf_ inner factor d= ddfm_ inner factor (\_ _ _ -> return 1) d

ddfm_ inner factor multF direction d= 
    Derivative_adj
        direction 
        (\_ _ _-> head $! multTerm factor $! inner d  ) 
        Center 
        (\st pos side -> multF st pos side d)
       
--drho_dt:: (Fractional a) => Reader (ValSet a) (Term a)        
drho_dt = d_ [Density] Time

--dp_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dx = d_ [Pressure] X

--dp_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dy = d_ [Pressure] Y

--dp_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dz = d_ [Pressure] Z

--drhou_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhou_dt =d_ [Density, U] Time   

--drhov_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhov_dt =d_ [Density, V] Time

--drhow_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhow_dt =d_ [Density, W] Time

drhodye_dt =d_ [Density, Dye] Time

--drhoT_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhoT_dt = df_ [Density, Temperature] specificHeatCv Time 

--dmewu_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dx = d_ [Mew, U] X

--dmewu_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dy = d_ [Mew, U] Y

--dmewu_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dz = d_ [Mew, U] Z

--dmewv_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dx = d_ [Mew, V] X

--dmeww_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dx = d_ [Mew, W] X

--dmewv_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dy = d_ [Mew,V] Y

--dmewv_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dz = d_ [Mew,V] Z

--dmeww_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dy = d_ [Mew,W] Y

--dmeww_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dz = d_ [Mew,W] Z

divergence vector = map ( uncurry dd_) $! zip vector (enumFrom X)  
    
gradient properties constantFactor = map (df_ properties constantFactor) $! enumFrom X  

integrateTerms integrate env =map (\term -> integrateOrder integrate Spatial Temporal [term env] )  

divergenceWithProps properties = divergenceWithPropsFactor properties 1  

divergenceWithPropsFactor properties constantFactor = 
    map (\x -> df_ (properties++[d2C x]) constantFactor x ) (enumFrom X)
 
divGrad properties constantFactor = divergence $ gradient properties constantFactor  

integrateOrder :: forall (m :: * -> *) b t.
                        Monad m =>
                        (t -> b -> m b) -> t -> t -> b -> m b
integrateOrder integrate i1 i2 term = integrate i1 term >>= integrate i2

multProps_adj props schemeType= 
    let func = prop_adj schemeType
    in foldl' 
        (\prev next pos side vs-> prev pos side vs *  func next pos side vs)          
        (\_ _ -> return 1)
        props  
     
squareDerivative properties constantFactor direction = 
    let foldedProps = multProps_adj properties
    in [ ddf_ (d_ (properties++properties) direction ) (constantFactor/2) direction 
        ,ddfm_ (d_ properties direction ) (constantFactor * (-1)) foldedProps direction ] 

pairedMultipliedDerivatives props1 props2 dir1 dir2  = 
    let p1 = multProps_adj props1 
        p2 = multProps_adj props2
    in [dd_ (d_ (props1++props2) dir1 ) dir2 ,
        ddfm_ (d_ props1 dir1 ) (-1) p2 dir2 ,
        ddfm_ (d_ props2 dir1 ) (-1) p1 dir2 ]

integUnknown :: AdjGraph
                  -> Property
                  -> DimensionType
                  -> [Term Double]
                  -> AdjPos
                  -> [Term Double]
integUnknown env unknownProp dimetype terms cellposition = 
    integ_adj env dimetype terms cellposition unknownProp


