# https://ferrite-fem.github.io/Ferrite.jl/stable/devdocs/FEValues/#devdocs-fevalues

struct SurfaceCellValues{T <: AbstractFloat, QR <: QuadratureRule{RefTriangle}, IP #, GIP
    } <: Ferrite.AbstractCellValues
    N_ref :: Matrix{T}  # scalar shape values φᵢ_hat(ξ_q)
    ∇N_ref:: Matrix{Vec{2,T}}   # reference gradients ∇φᵢ_hat(ξ_q)

    N_surf :: Matrix{T} # shape values on surface
    ∇N_surf:: Matrix{Vec{3,T}} # tangential surface gradients JG⁻¹∇φᵢ
    detJdV :: Vector{T} # √det(G)×w_q per quadrature pt.

    #M_geo :: Matrix{T}  # geometric shape values

    qr :: QR
    ip  :: IP
    #gip :: GIP
    mapping :: SurfaceTriangleMapping    
end

function SurfaceCellValues(
    ::Type{T},
    qr::QR,
    ip::LTI
) where { T <: AbstractFloat,
    QR <: QuadratureRule{RefTriangle}}
    #gip = ip
    no_q = length(Ferrite.getweights(qr))
    no_bfs = getnbasefunctions(ip)
    mapping = SurfaceTriangleMapping()

    # initial allocation for reference-space arrays
    N_ref = zeros(T, no_bfs, no_q)
    ∇N_ref = zeros(Vec{2,T}, no_bfs, no_q)

    # Compute reference values and gradients at every quad.pt
    for (q,ξ) in enumerate(Ferrite.getpoints(qr))
        for i in 1:no_bfs
            N_ref[i,q] = T(Ferrite.reference_shape_value(ip, ξ, i))
            ∇N_ref[i,q] = Vec{2,T}(Tuple(Ferrite.reference_shape_gradient(ip,ξ,i)))
        end        
    end
    #M_geo = copy(N_ref)

    #initalise surface-space arrays (filled by reinit!)
    N_surf = zeros(T, no_bfs, no_q)
    ∇N_surf = zeros(Vec{3,T}, no_bfs, no_q)
    detJdV = zeros(T, no_q)

    return SurfaceCellValues{T, QR, typeof(ip) #, typeof(gip)
    }(
        N_ref, ∇N_ref,
        N_surf, ∇N_surf, detJdV,
        #M_geo,
        qr, ip,# gip,
        mapping
    )

end

# T defaults to Float64
SurfaceCellValues(qr::QuadratureRule{RefTriangle}, ip::LTI) = SurfaceCellValues(Float64, qr, ip)

# Ferrite interface methods

Ferrite.getnquadpoints(cv::SurfaceCellValues) = length(cv.detJdV)
Ferrite.getnbasefunctions(cv::SurfaceCellValues) = size(cv.N_surf, 1)
#Ferrite.getngeobasefunctions(cv::SurfaceCellValues) = size(cv.M_geo, 1)

Ferrite.get_fun_values(cv::SurfaceCellValues) = cv.N_surf
#Ferrite.get_fun_gradients(cv::SurfaceCellValues) = cv.∇N_surf

#Ferrite.get_geometric_values(cv::SurfaceCellValues) = cv.M_geo
#Ferrite.geometric_value(cv::SurfaceCellValues, q::Int, i::Int) = cv.M_geo[i,q]

Ferrite.getdetJdV(cv::SurfaceCellValues) = cv.detJdV
Ferrite.getdetJdV(cv::SurfaceCellValues, q::Int) = cv.detJdV[q]

Ferrite.shape_value(cv::SurfaceCellValues, q::Int, i::Int) = cv.N_surf[i,q]
Ferrite.shape_gradient(cv::SurfaceCellValues, q::Int, i::Int) = cv.∇N_surf[i,q]

Ferrite.shape_value_type(cv::SurfaceCellValues{T}) where {T} = T
Ferrite.shape_gradient_type(cv::SurfaceCellValues{T}) where {T} = Vec{3,T}

Ferrite.function_interpolation(cv::SurfaceCellValues) = cv.ip
Ferrite.geometric_interpolation(cv::SurfaceCellValues) = LTI()
Ferrite.get_mapping(cv::SurfaceCellValues) = cv.mapping

function Ferrite.reinit!(
    cv :: SurfaceCellValues{T},
    cell:: Ferrite.CellIterator
) where {T}
    reinit!(cv, collect(Ferrite.getcoordinates(cell)))
    return cv
end

function Ferrite.reinit!(
    cv::SurfaceCellValues{T},
    cellnodes::AbstractVector{<:Vec{3}}
) where {T}
    nodes = collect(Vec{3,T}, cellnodes)
    #no_q = getnquadpoints(cv)
    no_bfs = getnbasefunctions(cv)
    for (q, (ξ , w)) in enumerate(zip(Ferrite.getpoints(cv.qr), Ferrite.getweights(cv.qr)))
        mv = compute_surface_mapping(nodes, ξ) #->J, G, invG, detJ
        cv.detJdV[q] = T(get_detJ(cv.mapping, mv)) * T(w)
        for i in 1:no_bfs
            cv.N_surf[i,q] = cv.N_ref[i,q]
            cv.∇N_surf[i,q] = map_shape_gradient(cv.mapping, cv.∇N_ref[i,q], mv)
        end
    end
    return cv
end


